# splitLayers: I1.3 — reentrant (concave) corner columns (2026-06-20)

**Track:** Phase 5 / WP 5.2 / **I1.3**, the *concave* half. Complements
`2026-06-17-i13-opposed-strip-ridge-consensus.md` (the *convex* ridge half) and builds on
`2026-06-16-shared-rail-consensus.md` (I1.1 S1/S2b, in code) +
`2026-06-16-s2-corner-bisector-consensus.md`.
**Status:** design, pre-implementation. **Needs David's approval before any C++**, and — per
the standing rule — the cell-closure must be **PROVEN in a sandbox model before the engine is
touched** (see D0). WSL-built (`zefra-foam`, OF v2512).
**Goal:** recover the bracket's concave-corner cluster (`skippedWarp 149` + net
`skippedConflict 189`, co-located at the reentrant corners) with full layers, checkMesh-clean —
the cluster that takes the bracket from 0.895 to ≥0.95.

## 1. What is PROVEN (do not re-litigate)

- **Per-rail warp strips (I1.4 W1) are non-viable in isolation.** Build-class instrument,
  2026-06-20: `openCellClass | warp 568 | recovered 0 | normal 0 | unattributed 0` — **every**
  open cell was owned by a relaxed-warp column. Quality was fine (skew 3.7, nonOrtho 59); the
  failure is pure **cell non-closure**. W1 reverted.
- **Single-consensus-assign is also non-viable** — forcing one assign per warp column revives the
  2026-06-12 assign-mismatch skew defect, which is the very reason the `(b)` warp skip exists.
- **Mechanism (grounded):** the warp/conflict columns are **co-located** at the reentrant corners.
  In the clean baseline the whole cluster skips together (clean). Building any **subset** leaves
  built-layer-stack ↔ skipped-original-cell (and divergent per-rail) interfaces that don't close
  (the W1 `OPPWALL colA(built)/colB(skipped)` samples). So the cluster must be recovered **as a
  coherent unit**, not column-by-column.

## 2. The principle (reused from the convex-ridge I1.3 + S2b)

Both prior wins close adjacent columns by making them **share ring labels** instead of each
inventing divergent rings: the convex-ridge doc shares **one consensus ring column at the ridge
wall point**; S2b shares the **canonical prefix ring labels** along a corner-bisector rail. The
concave corner is the same idea applied to a **cluster**: the columns meeting at a reentrant
corner must share a **consensus corner-rail fan** so their layer cells close against each other
and against the corner, rather than each column carving independently.

This is the plan's "reentrant corner wall cells get a **fanned column — shared corner rails,
layers wrapping the corner** (the standard commercial approach)."

## 3. Approach (to be PROVEN before coded)

1. **Cluster identification.** Detect reentrant-corner cells/points: wall cells with >1 patch
   face, and shared wall edges/points that are **concave** (the face-centre test already used,
   inverted from the convex test at splitLayers.C L850/L1158). Group columns sharing a reentrant
   corner point into one **cluster** (union-find over reentrant-corner adjacency — mirrors I1.1's
   rail-sharing graph, but keyed on the concave corner).
2. **Consensus corner-rail fan.** Per cluster, build **one** shared ring set at the reentrant
   corner (placed along the smoothed inward bisector at the corner point, as the convex case does
   along the ridge normal), and have every column in the cluster **reuse those ring labels** on
   their corner-side rails (registry keyed by corner point, like `ridgeRings`). Layers wrap the
   corner on the shared fan.
3. **Coherent build-or-terminate.** A cluster is built **only if the whole cluster closes**;
   columns that cannot join the fan cleanly **terminate on the shared corner ring** (not skip in
   isolation, which is what broke W1). No partial-subset builds.
4. **Strips/caps unchanged downstream** — once the cluster shares ring labels, the existing
   PASS-3/4 strip + cap construction consumes them (same as the convex-ridge fix); flat/non-corner
   columns are **byte-identical**.

## 4. D0 — PROVE-FIRST: the sandbox cell-closure model (BEFORE any CONSTRUCTION edit)

**This step is mandatory and comes before any CONSTRUCTION (mesh-changing) edit (C1).** It is
offline Python — no OpenFOAM, no engine edit at all — and is the direct guard against repeating
W1, which I authored without proving closure. (The C0 detector in §5 is behaviour-neutral
instrumentation — the S0/W0 category — so it is not a "modify by an untested theory"; even so, D0
is sequenced FIRST in §5 so the closure proof gates everything that follows.)

- **Inputs:** representative reentrant-corner geometry — synthetic (a parametric concave bracket
  corner) **and** the *real* bracket-corner columns, whose rail + ring coordinates we already have
  from the `OPPWALL|A/B|…|rail|…|seg0|k:(x y z)` dumps in the W1 log.
- **Model:** implement the §3 fan construction (shared corner-rail ring labels + wrap) as a pure
  geometric/topological model — points, rings, the per-cell faces (caps + strips), and the
  cell→face incidence.
- **Verdict:** for **every** layer cell, compute the signed-face-area sum `|Σ ±A_f| / max|A_f|`
  (the exact checkMesh open-cell criterion) and require it ≈ 0 for the whole cluster; also check
  no face is double-owned and winding is consistent. Sweep the corner angle (down to acute) and
  the per-rail length divergence (the warp severity) to find where closure holds.
- **Gate:** only when the model proves closure across the swept corner/warp range do we author the
  C++. If a regime won't close in the model, it stays a documented skip — decided on the model,
  not in the engine.

(This mirrors the stack-index sandbox check that correctly de-risked W1's indexing — that part
never crashed. The piece W1 *didn't* model was closure; D0 fixes that gap.)

## 5. Staging (each gated on a fixture; sweeps=0, edge placement)

- **D0 — the sandbox cell-closure model (§4). FIRST.** Offline; no engine edit. No construction is
  authored until it proves closure across the swept corner/warp range.
- **C0 — reentrant-cluster detector, emit-only** (behaviour-neutral instrument, S0/W0 category):
  `SPLITLAYERS|cornerCluster|<patch>|clusters|<n>|cells|<n>`; verify it lands on the bracket corner
  cluster (and is ~0 on flat externals / ahmed). May run alongside D0; required before C1.
- **C1 — consensus corner-rail fan + coherent build-or-terminate** (the construction D0
  validated). Gate: bracket corner cluster closes (open cells 0 via the build-class instrument),
  warp+conflict residual drops, checkMesh-clean, **ahmed/cow/wing byte-identical** (no reentrant
  clusters there beyond their own corners → guard with the no-regression check).
- **C2 — compose with the convex-ridge fan** (cow ears) and S2b so all three corner mechanisms
  coexist; re-run the suite census.

## 6. Verification / DoD

- Bracket: corner-cluster open cells → 0; `skippedWarp` + net `skippedConflict` fall enough to
  clear **≥0.95 thickness coverage**; checkMesh-clean (or orchestrator-clean for genuine residual
  skew only); nonOrtho ≤ 70.
- lbend family + ahmed/cow/wing: **no regression** (byte-identical where no reentrant cluster).
- The build-class `openCellClass` instrument stays on as the guard; `report_fingerprint` stable (I4.3).

## 7. Risks & the standing discipline

- **PROVE BEFORE MODIFY ([[feedback_prove_before_modify]]).** No CONSTRUCTION (mesh-changing) edit
  until D0's model proves closure (C0 is a behaviour-neutral instrument, sequenced after D0 anyway).
  State every cause as a hypothesis until measured. (This doc exists because W1 was built on an
  unproven closure assumption.)
- **Touches ring creation (PASS 3)** — the core. Mitigated: corner clusters are a small set;
  flat/interior placement unchanged (byte-identical externals); the open-cell + build-class
  instruments catch any new non-closure at the source.
- **Concave vs convex discrimination** must be exact (reuse the face-centre test; convex ridges go
  to the 2026-06-17 doc, concave corners here).
- **Determinism:** union-find + corner-keyed registry filled in column order, tie-broken by lowest
  point/column id → order-independent; assert via the report fingerprint.
- **Coherent-or-terminate, never partial-subset** — the W1 lesson: a built column abutting a
  skipped cluster neighbour does not close. The cluster builds whole or its unjoinable members
  terminate on the shared corner ring.

## 8. Open decisions (for review)

- **D1 — cluster scope:** group by shared reentrant-corner POINT (recommended; mirrors the convex
  ridge-point keying) vs by corner CELL (>1 patch face). Point-keying composes with the convex fan.
- **D2 — terminate vs skip for unjoinable members:** terminate on the shared corner ring
  (recommended — keeps coverage, closes cells) vs leave skipped (coverage-only fallback if the
  model shows termination skews).
- **D3 — build the D0 model on synthetic corners only, or also feed the real OPPWALL-dumped
  bracket coordinates** (recommended: both — synthetic for the sweep, real for ground truth).
- **D4 — sequence vs I2.1:** I1.3 is the bracket's ≥0.95 path but is the hardest code; if the D0
  sweep shows acute corners won't close, bank the bracket at its documented-skip floor and let
  I2.1 (non-hex: ahmed/cow/wing) carry near-term suite coverage. Decide on the model's result.
