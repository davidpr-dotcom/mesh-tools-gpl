# splitLayers: S2 — corner-bisector (joints) conflict consensus (I1.1 design, 2026-06-16)

**Track:** Phase 5 / WP 5.2 / I1.1, **S2** — the last conflict lever after S1.
**Reads with:** `2026-06-16-shared-rail-consensus.md` (the S0/S1/S2/S3 frame).
**Status:** design + a diagnostic FIRST (the S0 discipline), pre-implementation.
This is the hardest geometry code in the track; it overlaps **I1.3 corner columns**.

## State after S1 (WSL 2026-06-16, convex skip default on)

bracket: cols 1593, **skippedConflict 202**, skippedWarp 73, skippedConvexRidge 238
(orchestrator-handled), checkMesh OK. The 202 is **joints-dominated** (S0: joints
192) — the corner-bisector case S1 cannot touch (S1 only harmonizes scale/assign
on *genuinely-shared* rails; these are *different* rails). ahmed/cow/manifold are
essentially clean post-S1 (+orchestrator). So S2 ≈ the bracket coverage ceiling.

## The mechanism (grounded)

`railsOf` is `EdgeMap<RailData>` keyed by `railKey = edge(rails[c][0], rails[c][1])`
(wall point + first interior joint). At a reentrant/convex **corner bisector**, two
adjacent columns A,B grow a rail from the **shared corner wall point** and reach the
**same first interior joint** → **same key** — but then their rails go through
**different cells** → **different full joint paths** (`rails[c]` differ from seg1
onward). `railsOf` can hold only ONE `RailData` per key (it stores A's full path +
rings). B's conflict check finds `it.val().joints != B.rails[c]` → **skip B**
(`reason="conflict"`). Geometrically: A and B genuinely **share segment 0** (the
first edge) but **diverge after it**; the per-whole-rail model can't represent
"shared prefix, divergent suffix", so B is dropped.

(Note: this is distinct from the PASS-3 `edgeClaim`/`nEdgeConflict` late skip —
that catches rails from *different* wall points merging onto one edge. The
corner-bisector is caught earlier, in the PASS-2 full-path conflict check.)

## S2-S0: diagnostic FIRST (instrument before designing the fix)

Before choosing the fix we must know the *structure* of the 192 joints-conflicts —
a single global angle/threshold won't have told us, just as it didn't for the
convex ridge. Add a behaviour-neutral instrument to the conflict check: when a
column is skipped on the **joints** arm, emit the **shared-prefix depth** (how many
leading rail points A and B share before diverging) and the **corner type**
(reentrant vs convex, from the wall-edge normals). Re-run the bracket.

Read-out decides the fix:
- **mostly share only seg0 (prefix depth 1), then diverge** → the **per-segment
  consensus** below recovers them cleanly (they only need to agree on seg0).
- **deeper shared prefixes / true multi-face corners** → the **corner-column fan**
  (I1.3) is required; per-segment consensus handles the shallow ones and the rest
  is I1.3.
- **reentrant vs convex split** → tells us whether the wrap is concave-bisector or
  convex-ridge-like (the convex ones may already be handled once the orchestrator
  disables the geometric ridge skip).

Patch sketch (mirrors the S0 instrument): in the PASS-2 conflict block, when
`it.val().joints != col.rails[c]` is the (or a) cause, compute the common-prefix
length of the two `labelList`s and emit
`SPLITLAYERS|jointsConflict|<patch>|prefix|<k>|cornerType|<reentrant|convex>`.

## Diagnostic RESULT (bracket, WSL 2026-06-16)

`SPLITLAYERS|jointsConflict|prefixHistogram|2:125|4:70|5:4` (199 ≈ the 202
conflict). **All are shared-prefix → divergent-suffix**: 125 share only seg0
(depth 2), 74 share a 3–4-segment prefix then diverge (depth 4–5). **None are
deep tangled multi-face corners → the I1.3 fan is NOT needed for the bracket;
per-segment consensus recovers the whole population** (the deep ones just share
more segments before diverging). Projected bracket: 1593 → ~1792 cols (≈85%) from
S2b; + convex ridges via the orchestrator + warp (I1.4) → the ~95% target.
**Refinement the depth-4/5 cases force:** they diverge *after* sharing 3–4
segments, and S1 only harmonized scale on *identical-full-path* rails — so S2b
must also harmonize scale on the **shared prefix segments** (S2b = S1's min-scale
generalized from whole-rail to per-shared-segment-group).

## Candidate fix — per-segment ring consensus (CONFIRMED approach)

Rings already live per **segment edge** in `segRings` (`EdgeMap<labelList>`); the
*conflict* and *creation*, though, are keyed per **whole rail** (`railsOf`). The
fix decouples them:

1. Negotiate rings **per shared segment**, not per whole rail. Two columns sharing
   segment 0 share *one* seg-0 ring set (consensus scale/assign on that edge, the
   S1 machinery applied per-segment); each then builds its **own** deeper segments
   on its divergent path.
2. Relax the PASS-2 conflict check from **full-path identity** to **shared-segment
   agreement**: a column is compatible if, on every segment it *shares* with an
   already-recorded rail, the ring set matches; divergence *after* the shared
   prefix is allowed (it's a different segment, its own rings).
3. Keep the strict skip as the fallback for the genuinely-irreconcilable (true
   multi-face corners → I1.3 fan), so S2 is — like S1 — strictly coverage-adding,
   never a regression; the I4.2 matrix + the quality-retry orchestrator backstop
   any recovered column that turns skew.

This is the doc's "canonical joint-path / per-segment" S2. The exact data-structure
change (segment-keyed RailData, or a shared-segment ring registry both rails
consult) is finalized **after** the diagnostic confirms the prefix-depth structure.

## Staging & relationship to I1.3 / I1.4

- **S2a** — the diagnostic (above), bracket re-run.
- **S2b** — per-segment consensus for the shallow (seg0-only) corner bisectors —
  the bulk if the diagnostic shows prefix depth ≈ 1.
- **Residual** — deep/multi-face corners → **I1.3 corner columns** (fanned/wrapped
  shared rails, the standard commercial approach); warp skips → **I1.4**.
- Each stage gated on the **bracket_corner** fixture + checkMesh-clean; the
  orchestrator demotes any recovered column that turns skew (S1's manifold lesson).

## S2b implementation architecture (detailed, 2026-06-16)

**Key code finding that makes this tractable:** the **topology already consumes
rings per segment edge** — PASS 4's `withForeignRings` and the lateral pass walk
`segRings` (`EdgeMap<labelList>`, keyed by segment edge). The whole-rail constraint
lives only in two places: (b1) the PASS-2 conflict check (whole-rail identity), and
(b2) PASS-3 ring *creation* (`rd.rings` per `railKey`, one rail per first edge, in
`railsOf`). So we do NOT rewrite the topology; we change creation + the check.

**The model:** rings are identified **per segment edge**; a column's rail is a
sequence of segments and its `rd.rings[k]` reference whichever point sits on the
segment that layer k falls in. Two rails sharing a prefix share those segments'
ring **point labels** (not duplicates) — so the mesh is consistent — and each
builds its own divergent-suffix segments.

**Changes, pass by pass:**

1. **Scale (extend S1).** Today S1 unites only *identical-full-path* rails. Extend
   the union to rails that **share the first edge** (`railKey`) so prefix-sharers
   get one component scale (MIN). Shared-prefix segments then have identical arc +
   scale → identical layer positions → the shared ring points coincide and can be
   reused. (Components stay local — a corner point's ~4 columns — not a global
   cascade.)

2. **PASS 3 — per-segment ring creation (the core change).** Introduce a segment
   registry `EdgeMap<labelList> ringPtsOf` (segment edge → ordered ring labels,
   wall-side→interior; both rails traverse a segment the same way since rails grow
   outward). Per column, per rail: walk segments accumulating arc; group layers k
   (at `h_k = scale·cumHeights[k]`) by the segment they fall in; for each segment
   edge — **reuse** its labels from `ringPtsOf` if present, **else create** the
   points + register. Assemble the column's own `rd.rings[c]` (now **per-column**,
   not per-railKey). Register `segRings` exactly as today (topology unchanged).

3. **Validity guard (robustness — never silently wrong).** Before reusing a
   segment's labels, verify the reusing rail agrees on that segment: same scale,
   same cumulative-arc-to-segment, same layer set. The diagnostic says all sharing
   is *prefix* sharing (same arc) so this always holds for the bracket; the guard
   makes it robust to any future non-prefix sharing — on mismatch, **fall back to
   the strict skip** (coverage-only, never a bad mesh).

4. **PASS 2 — relax the conflict check.** Replace whole-rail `joints/assign`
   identity with: a column is compatible if, on every segment it *shares* with an
   already-created rail, the guard (3) passes; divergence after the shared prefix
   is allowed. Genuinely-irreconcilable → strict skip (fallback). nSeg is no longer
   a blocker (divergent rails legitimately differ in length).

5. **PASS 4 — unchanged.** It already reads `segRings` (per segment) for side
   faces and the per-column stack for cells; per-column `rd.rings` feeds the layer
   interfaces as before. Audit that nothing reads `railsOf[key].rings` assuming
   one-rail-per-key (the optimize pass `forAllConstIters(railsOf,…)` must iterate
   the new per-column store instead).

**Why this is the most-robust long-term choice (not a patch):** it removes the
incorrect whole-rail identity assumption at its root (rings belong to segments),
makes prefix-sharing fall out automatically with no special-casing, keeps the
strict skip as a coverage-only safety net, and leaves the topology untouched. It
generalizes S1 (whole-rail scale consensus → per-segment ring consensus) rather
than bolting on a corner special case.

**Staging:** the per-segment model handles all prefix depths uniformly, so it's one
coherent change rather than depth-2-then-depth-4. Build it behind the strict-skip
fallback so each WSL iteration is safe; gate on bracket_corner + the orchestrator.
True multi-face corners (if any surface) and warp remain I1.3 / I1.4.

## DoD

bracket `skippedConflict` materially down from 202 (target the joints share the
diagnostic says is seg0-recoverable; the plan's "< 50" is the combined S2+I1.3
target), checkMesh-clean, no regression on the I4.2 matrix, deterministic.

## Risks

- **Per-segment decoupling touches the ring/conflict core** — higher regression
  risk than S1. Mitigated by: the strict-skip fallback (coverage-only), the
  bracket_corner fixture, and the orchestrator quality backstop.
- **Determinism** — segment consensus must be order-independent (tie-break by
  lowest column id), asserted by the report fingerprint.
- **Overlap with I1.3** — don't over-build S2; recover the shallow bisectors, hand
  true corners to I1.3. The diagnostic draws that line.
- **GPL wall** — all executor mechanics, behind the D7 boundary, as before.
