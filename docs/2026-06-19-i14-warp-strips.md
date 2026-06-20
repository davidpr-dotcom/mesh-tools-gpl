# splitLayers: per-rail warp strips (I1.4 design, 2026-06-19)

**Track:** Phase 5 / WP 5.2 / **I1.4** — the measured top bracket lever after I1.1 (S1+S2b).
**File:** `mesh-tools-gpl/applications/splitLayers/splitLayers.C` (2278-line build).
**Goal:** let a column whose four rails **disagree on the ring→segment assignment** still build
its layers, by constructing the lateral strips and ring caps **per rail** with sub-faces that
**terminate on the joint** where a rail lacks a ring — instead of skipping the column. Apply the
same relaxation to the **S2b recovery warp-guard** so corner-bisector sharers blocked only by warp
are recovered too.
**DoD:** bracket warp skips < 30 (from the assignment-warp share of 149); net conflict falls as the
S2b warp-aborts (145) succeed; perturbed L-bend coverage ≥ 98%; checkMesh-clean; no I4.2 regression.
**Status:** design — for review BEFORE any C++. WSL-built (`zefra-foam`, OF v2512). Diagnostic-first
(W0 instrument) before construction code.
**Reads with:** `2026-06-16-shared-rail-consensus.md` (I1.1 S0→S2b, the seeding architecture this
extends), `Zefra/2026-06-12-splitlayers-industrialization-plan.md` (I1.4 track + DoD),
`Zefra/2026-06-17-phase5-implementation-sequence.md` (Stage B, diagnostic-first).

---

## 1. The measured motivation (2026-06-19 bracket diagnostic)

The bracket suite-eval, post-I1.1, with the new recovery census surfaced:

```
concave_bracket  conflict 320→net 189 (recovered 131)  warp 149
  recoverAbort {noCanon 0, warp 145, scale 44, edgeClaim 0}  ridgePts 314
```

Reading it: I1.1 (S1 min-scale + S2b corner-bisector recovery) already recovers 131 conflicts; the
conflict residual is **189**, and there are **no genuine corner orphans** (`noCanon 0, edgeClaim 0`).
The dominant blocker is **warp**:

- **149 standalone `skippedWarp`** columns (the `(b)` rail-assignment skip), and
- **145 S2b recovery aborts** with `why == "warp"` — conflict columns that *would* recover but their
  four rails disagree on assignment.

That is **up to ~294 faces gated on one mechanism** — by far the largest single lever, and enough to
clear the ≥95% bracket gate. (`ridgePts 314` is the convex-ridge / cow-ear class — I1.3, separate.)

**Important unknown → W0 first.** The standalone count `149` lumps **two** distinct causes that both
set `reason="warp"`: the `(a2)` **degenerate** chain (duplicate point among the four rails, L908–931)
and the `(b)` **assignment-disagreement** (L933–944). I1.4 fixes only `(b)`; `(a2)` is a genuine
degeneracy that should stay skipped. The S2b-abort `145` is pure `(b)` (it runs only the assignment
test, L1340–1347). So the addressable standalone share is unknown until we split the counter — hence
**W0 (a behaviour-neutral diagnostic) before any construction**, mirroring the S0 discipline.

## 2. The current mechanism (what enforces "all four rails agree")

Three places assume a column has **one** assignment table shared by all four rails:

1. **The `(b)` warp skip (L933–944).** `a0 = computeAssign(col, 0)`; for `c = 1..3`, if
   `computeAssign(col, c) != a0` → `skip`, `reason="warp"`, `++nWarp`. `computeAssign` (L777–797)
   walks each rail's segment arc-lengths and returns, per layer `k`, the segment index the layer
   height `scale·cumHeights[k]` lands in. Rails whose segments have **different arc-lengths**
   (a warped cell — non-parallel opposite edges) put layer `k` on segment `i` for one rail and
   `i±1` for another → the tables differ → skip.
2. **The ring-cap builder (PASS 4 §2, L1658–1670).** Each cap face at layer `k` is the quad
   `{ ringsFor(ci,col,3-c)[k] for c in 0..3 }` — it indexes **the same `k`** on all four rails,
   valid only if every rail has a ring at layer `k` (i.e. they agree).
3. **The lateral strip builder (PASS 4 §5, L1746–1913).** `cutK` = layers with `assign[k]==i` from
   the column's **single** `assign` (`assignFor`, L1559–1563 → `railsOf[railKey(col,0)].assign`).
   For each sub-cut `s`, `pickPt(c, lower)` returns `ringsFor(ci,col,c)[cutK[idx]]` on **both**
   bounding rails `ca, cb` — again the same layer index on both rails.

So the `(b)` guard is not cosmetic: it protects §2 and §5 from indexing a ring a rail does not have
on this segment (which would produce a twisted/degenerate face and, historically, the
`changeMesh` boundary-face crash the OOB/BADFACE instruments hunt). Relaxing the guard therefore
**requires** generalizing §2 and §5 to per-rail assignment first.

## 3. Root cause

`scale` is per-column but the **assignment is per-rail geometry**: two rails of the same column with
different segment arc-lengths legitimately distribute the same layer heights onto different segments.
The construction was written for the common case (all rails congruent → one table) and **skips** the
warped case rather than building a strip whose two sides step at different layers. The fix is a strip
that is the union of both rails' steps, **closed on the joint** wherever only one side has a ring.

## 4. The design — per-rail assignment + joint-terminated sub-faces

Keep `scale` per-column (consensus already set it); make **assignment per-rail**. Per rail `c`,
`assign_c = computeAssign(col, c)` (already computed in the guard — retain all four instead of
comparing-then-discarding). Then:

### 4.1 Lateral strips (§5 generalization — the 2026-06-12 "sub-face cuts terminate on joints")
For a side face bounded by rails `ca, cb`, the sub-faces currently step at the shared `cutK`. Replace
with a **merge walk** over the two rails' ring heights on this segment `i`:

- Gather `ringsCutA = { (h_k, ringsFor(ca)[k]) : assign_ca[k]==i }` and likewise `ringsCutB`.
- Walk both in height order. Each sub-face spans from the previous boundary to the next, where each
  side independently advances to its next ring **if it has one at this step**, else **stays on the
  joint** `col.rails[c][i]` / `[i+1]` (the lower/upper segment endpoint). The existing `pickPt`
  already terminates on the joint at `s==0` and `s==nSub-1`; this extends that to any step where a
  side lacks a ring — so a strip with `nA` cuts on one side and `nB` on the other emits
  `max(nA,nB)+1` sub-faces, the short side repeating its joint (a triangle there, not a twisted quad).
- The owner/neighbour stack indexing (`ringsBefore[i]+s`, the `chk`/`term` clamps, L1842–1871) keys
  off the **owning rail's** ring count on the segment; with per-rail counts, `ringsBefore` becomes
  per-rail (`ringsBefore_c`) and `s` maps through the side's own ladder. The `term` clamp (variable-n
  termination, L1588–1591) already handles a side running shorter than its stack — joint-termination
  is the same idea applied laterally.

### 4.2 Ring caps (§2 generalization)
A cap at layer `k` uses only the rails that have a ring at `k` on their segment; a rail that places
layer `k` on a different segment contributes its **joint** (segment endpoint) instead of a ring, so
the cap degenerates from a quad to a triangle/pentagon that still closes. Same `withForeignRings`
insertion applies. (This is the subtler half — caps are the cells' tops; a mis-built cap is the
classic open cell. Build behind `saneCheck` and the `self_check` refuse-guard.)

### 4.3 Apply to BOTH skip sites
- **Standalone `(b)` (L933–944):** for assignment-disagreement (not `(a2)` degeneracy), do **not**
  skip — let PASS 3/4 build per-rail. Keep `(a2)` degenerate skip unchanged.
- **S2b warp-guard (L1340–1347):** the recovery aborts when `computeAssign(col,c)!=a0`; with per-rail
  strips this guard is no longer needed for the assignment case — the recovered sharer carries its
  per-rail assigns through the same builder. (Keep the `scale`/`edgeClaim` guards: those are real.)

## 5. Staged rollout (each gated on a fixture; sweeps=0, edge placement)

0. **W0 — diagnostic FIRST (behaviour-neutral).** Split `nWarp` into `nWarpDegenerate` (a2) and
   `nWarpAssign` (b), and emit `SPLITLAYERS|warpCause|<patch>|degenerate|<i>|assign|<i>`. Re-run the
   bracket. This sizes the *addressable* standalone share (only `assign`) and sets the realistic DoD
   (warp < 30 is meaningful only against `nWarpAssign`, not the degenerate floor). Mirrors S0; the
   parser already ignores unknown `SPLITLAYERS|` records, and the meshcore report parser will surface
   `warpCause` with the same additive treatment used for the recovery census.
1. **W1 — per-rail lateral strips** (§4.1) behind the standalone `(b)` skip, *plan-then-commit*
   (compute all sub-faces + run `saneCheck` with **zero** mesh side effects, then emit) exactly like
   S2b. Gate: bracket `nWarpAssign` drops sharply, checkMesh-clean, no quality regression, OOB/BADFACE
   counts stay 0.
2. **W2 — per-rail ring caps** (§4.2). Gate: bracket warp < 30, net conflict unchanged or better,
   perturbed L-bend ≥ 98%, checkMesh-clean.
3. **W3 — relax the S2b warp-guard** (§4.3) to route recovered sharers through W1/W2. Gate: the 145
   warp-aborts recover (census `recoverAbort.warp` → ~0), net conflict drops accordingly, I4.2 green.

Each stage is independently shippable and measured (the suite-eval census is the scoreboard).

## 6. Decisions for review (recommendations)

- **D1 — joint-termination vs ring-insertion at the short side.** *Recommend* joint-termination
  (the short side rests on the segment endpoint), per the 2026-06-12 sketch: it never invents a ring
  the rail's geometry doesn't support, so y⁺ stays honest. Inserting an interpolated ring on the short
  side would equalize counts but fabricate height — rejected.
- **D2 — keep the all-agree fast path.** *Recommend* yes: if `assign_c == a0` for all `c` (the common
  case, every congruent cell), use today's exact code path → byte-identical output for all current
  clean columns (determinism + zero regression on externals). Per-rail is the *else* branch only.
- **D3 — W2 ring caps scope.** *Recommend* gate W2 on cow_ear staying clean and the bracket; if a cap
  degeneracy proves brittle on the tightest warp, fall back to leaving those specific columns skipped
  (coverage-only) rather than shipping a skew cap — the orchestrator + `self_check` enforce that floor.
- **D4 — DoD denominator.** Use `nWarpAssign` from W0 as the < 30 target denominator; report the
  degenerate floor separately so a residual `(a2)` count isn't read as an I1.4 miss.

## 7. Risks & past mistakes carried in (§6 lineage)

- **Transactional construction only.** Plan every sub-face + cap and run all `saneCheck`/bounds guards
  with **zero** side effects, then commit `addFace`/`modifyFace` — a pre-abort `addPoint`/`addFace`
  orphans points → "unused points"/dangling-face checkMesh fail (the S2b model, L1353–1413).
- **Keep the instruments on.** `saneCheck` (size-neighbourhood), `chk`/`term` (stack OOB), `badFace`
  (dangling owner / `nei<0 && patch<0`), `latRemap` (never-crash) must stay green; per-rail indexing
  is exactly where an off-by-one re-introduces the `changeMesh` boundary-face crash.
- **Degenerate ≠ warp.** Do not let W1 try to build the `(a2)` duplicate-point columns; that is a real
  degeneracy (converging rails), not an assignment mismatch. W0 separates them.
- **Determinism (I4.3).** The merge walk must be order-independent (walk by height then by lowest
  point label); assert the report fingerprint is byte-identical across two runs.
- **A correlate is not a discriminator.** Measure with checkMesh + the I4.2 comparator on the BUILT
  mesh, not by eye; the saneCheck spread is a debug aid, not the quality gate.
- **GPL wall.** All of this is executor mechanics inside GPL `splitLayers`, behind the D7 process
  boundary; no planner intelligence moves here.

## 8. Fixtures (capture into the I4.2 matrix)

- **bracket_corner** (already captured): must show `nWarpAssign` drop W1→W2 and the S2b warp-aborts
  recover at W3; checkMesh-clean, nonOrtho ≤ 70 throughout.
- **perturbed L-bend** (warp stress; per the I1.4 plan DoD): coverage ≥ 98% after W2.
- **cow_ear + lbend family + cavity:** stay green (no regression) at every stage.

## 9. Build + run (David, WSL — `zefra-foam`, OF v2512) — once each stage's code lands

```bash
cd /mnt/c/Users/traca/Desktop/fullCrossflow/mesh-tools-gpl/applications/splitLayers && wmake
cd /mnt/c/Users/traca/Desktop/fullCrossflow/Zefra
python3 tools/2026-06-18-layer-suite-eval.py concave_bracket    # census line is the scoreboard
# W0: grep the new warpCause split
grep "SPLITLAYERS|warpCause" myCases/wp54_spike/concave_bracket/layer_suite/log.splitLayers
python3 -m pytest meshcore/tests -q                              # parser/report contract tests
```
(Exact per-stage commands restated when each patch is handed over.)

## 10. Files to be changed (none overwritten without showing first)

- `mesh-tools-gpl/applications/splitLayers/splitLayers.C` — W0 counters + `warpCause` emit; then the
  per-rail strip (§4.1), per-rail caps (§4.2), and S2b guard relaxation (§4.3), staged W1→W3.
- `Zefra/meshcore/executors/splitlayers_report.py` — additive `warpCause` parse (same pattern as the
  2026-06-19 recovery-census change) + a test.
- `mesh-tools-gpl/tests/regress/` — perturbed-L-bend warp fixture for the I4.2 matrix.

**First concrete step on approval:** I author the **W0 warp-cause instrument** (behaviour-neutral,
like S0) for you to `wmake`/run, so the degenerate-vs-assignment split sets the real W1 target before
any construction code is touched.
