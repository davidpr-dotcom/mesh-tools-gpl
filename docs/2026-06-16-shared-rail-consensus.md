# splitLayers: shared-rail consensus (I1.1 design, 2026-06-16)

**Track:** Phase 5 / WP 5.2 / **I1.1** — the single biggest coverage lever.
**Goal:** replace the *first-come-wins* rail recording with a **consensus pass** so
neighbouring columns that share a rail agree by construction and **keep their
layers** instead of being skipped — and so two columns meeting on a sharp **convex
wall edge** (cow ears) build the seg0 strip **once, correctly oriented**, retiring
the interim `skippedConvexRidge`.
**DoD:** bracket `skippedConflict` < 50 (from 345); `skippedConvexRidge` → 0 with
the cow ears split **checkMesh-clean**; no quality regression on the I4.2 matrix.
**Status:** design — for review BEFORE C++. WSL-built (`zefra-foam`, OF v2512).

---

## 1. The current mechanism (what we are replacing)

Lines below are `applications/splitLayers/splitLayers.C` (1401-line build).

**Data model.** A `Column` (L159) is a wall face + a `chain` of cells + four
`rails` (`FixedList<labelList,4>`, one point-path per quad wall vertex, length
`m+1` for a chain of `m` cells) + a `scale ∈ (0,1]`. A `RailData` (L173) is a
rail's record: its `joints` (the **FULL point path = identity**), `rings` (new
point labels, filled PASS 3), `assign` (segment index per ring), `nSeg`, `scale`,
`si`.

**PASS 1 — identification (L319-418).** Per wall face, claim the owner cell
(`claimed[own] = columns.size()`, **first-come**), grow the four rails through the
cell via `railStep`, and extend the chain through internal faces while
`minRailLen() < needed[si]` and the next cell is unclaimed/hex/non-patch. Then
`col.scale = min(1, minRailLen()/needed[si])` (L415) — a short column (its onward
cell was claimed or non-hex) gets a **smaller scale** than its neighbour.

**PASS 2 — conflict pre-passes (L420-695).** `railsOf` is an
`EdgeMap<RailData>` keyed by `railKey` = the **first segment edge**
`edge(rails[c][0], rails[c][1])` (L426). Iterating columns **in order**, each is
checked and, if it survives, **records** its rails:

- **conflict (a), L513-532:** for each rail, if `railsOf[railKey]` already exists
  and disagrees on `nSeg`, `scale` (>1e-6 rel), `assign`, **or** the full `joints`
  path → `skip`/`"conflict"`, `++nConflict`.
- degenerate (a2) / warp (b) / cross (c) — other strict skips.
- **record (L681-694):** `if (!railsOf.found(railKey)) railsOf.insert(key, rd)`
  — **the first column to reach here fixes the rail's `(scale, assign, joints,
  nSeg)`; every later column sharing that key validates against it.**

**convex-ridge skip (L464-505, I0.1 interim).** Before (a): if a wall-face edge is
sharply convex (`convexAngle`, default 70°) and the neighbour cell is claimed by a
different column → `skip`/`"convexRidge"`, `++nConvex`, for **both** columns. This
is the placeholder I1.1 retires.

## 2. Root cause (restated)

First-come-wins is *order-dependent strictness*: whichever column is visited first
imposes its `(scale, assign, joints)` on a shared rail; a neighbour that legitimately
diverges is **skipped** rather than reconciled. That is the bulk of `skippedConflict`
(bracket 345, manifold 233k).

**Conflict taxonomy (sharpens what each stage targets).** `scale` is a **per-column**
property — `min(1, minRailLen/needed)` over the column's *four* rails (L415) — so two
columns that share **one** rail (identical `joints`, hence identical `nSeg`) can still
carry **different `col.scale`**, because their *other* rails differ in length. On the
shared rail the conflict check then trips on `scale` alone. **This is the dominant,
cleanly-recoverable case → S1 (min-scale consensus).** A second, harder case: two
columns share only the **first edge** but diverge afterward (corner bisectors) —
genuinely *different* rails colliding on the `EdgeMap` first-edge key, so only one
`RailData` can be stored and the other conflicts (`joints != col.rails[c]`). That is
**not** a scale problem; it needs key disambiguation / canonical-path handling → **S2**
(and the residual goes to I1.3 corner columns). The **convex ridge** is the S2/S3
family forced by a sharp convex wall edge: the two seg0 strips meet as `strips:opposed`
and come out wrong-oriented/skew (I0.1) — skipped today until S3's shared seam lands.

## 3. The consensus design

Insert a **negotiation phase between identification and validation**. Instead of
the first column dictating, gather all rail instances, group the ones that *should*
agree, negotiate one consistent record per group, seed `railsOf` with it, then let
the existing checks validate against the **negotiated** record (so agreement is the
norm and a skip is the genuine exception).

### 3.1 The rail-sharing graph

- **Nodes:** rail instances `(ci, c)` — column `ci`, rail index `c ∈ {0..3}`.
- **Edges (same-key):** two instances share an edge if they have the same
  `railKey` (first segment edge) — today's collision, the conflict source.
- **Edges (convex-ridge seam):** two instances are **also** joined if their wall
  vertices sit on opposite sides of a **sharp convex wall edge** (the L464 test),
  i.e. the cow-ear case where the keys *differ* (distinct crease points) but the
  strips must share one negotiated seam. This is the new adjacency that lets the
  ridge be built once.
- **Components:** union-find over those edges → a *rail group* is one connected
  component that must share a record.

### 3.2 Per-component negotiation

For each component, compute ONE record:

1. **scale = min over the component** (the most-constrained column sets the
   achievable height — every member uses it, so no scale-mismatch conflict), then
   **Lipschitz-smooth** along the wall so neighbouring groups' scales never jump by
   more than the grading bound (the spec `expansionRatio`-derived growth; cf.
   `meshcore.sizing.grading.check_lipschitz`). Min-then-smooth keeps the field
   monotone-safe and deterministic.
2. **canonical joint-path:** choose one `joints` per shared rail — the path of the
   member with the **most cells** (the longest/most-constrained chain, a superset
   that the shorter members' rings can map onto). Shorter members re-segment
   against this canonical path under the shared scale.
3. **assign / nSeg:** recompute `assign` from the shared `scale` + canonical path
   (reuse `computeAssign`); `nSeg` reconciled per §6 decision (recommended:
   **min-nSeg** on the shared rail — both members layer to the shallower common
   depth at the seam, keeping layers over depth).

Seed `railsOf` with the negotiated record for every key in the component **before**
the validation loop. The existing (a) check then passes for all members (they match
the consensus), so they keep their layers.

### 3.3 Convex-ridge seam

Because the ridge adjacency (§3.1) put the two ear columns in one component, the
negotiation produces **one shared seg0 (and deeper) ring path along the crease**,
oriented consistently (outward of the convex edge). Both columns build their seg0
strip against that **single** seam instead of each building an opposed strip — the
`strips:opposed` construct never forms. If the simple shared seam is not enough for
the tightest ear curvature to go checkMesh-clean, the residual is the **fanned/
wrapped corner column** of **I1.3** (which consumes the same shared-rail record);
until I1.3, the `skippedConvexRidge` fallback stays as the safety net (§5).

### 3.4 Re-validate + never regress

After seeding, run the existing conflict/degenerate/warp/cross checks unchanged.
Any column the consensus could **not** reconcile (e.g. a genuinely incompatible
joint topology) still hits the existing skip — so I1.1 is **strictly a recovery**:
it can only turn former skips into kept layers, never the reverse. The I4.2 matrix
(bracket_corner, cow_ear, lbend family) gates that no fixture regresses.

## 4. Where it hooks (implementation shape)

Two-phase PASS 2:

```
PASS 2a  GATHER     build rail-instance list + railKey index + convex-ridge pairs
         GRAPH      union-find components over same-key ∪ ridge-seam edges
         NEGOTIATE  per component: min-scale → Lipschitz smooth → canonical joints
                    → assign/nSeg; write negotiated RailData into railsOf (seeded)
PASS 2b  VALIDATE   the existing loop (conflict/degenerate/warp/cross), now
                    comparing against the SEEDED consensus; record is already
                    present so first-come no longer matters; residuals skip as today
```

The change is localized: `railsOf` is **pre-populated** by 2a; the L681-694
"record if absent" becomes a no-op for negotiated keys (or is removed). `Column`
gains nothing; a small `negotiatedScale[ci]` / per-rail override may be needed if a
column's own `scale` must be lowered to the component min before `computeAssign`.

## 5. Staged rollout (each stage gated on a fixture, sweeps=0, edge placement)

0. **S0 — diagnostic FIRST (instrument, don't dial blind).** Before any consensus
   code, split the single `conflict` counter into its **sub-causes** (`nSeg`,
   `scale`, `assign`, `joints` mismatch) and re-run the bracket + manifold. This
   tells us whether the 345 (and 233k) skips are dominated by **same-rail/
   different-scale** (→ S1 is the lever) or **different-joints/same-first-edge**
   corner bisectors (→ S2 keying/corner is the lever). Patch:
   `2026-06-16-i11-s0-conflict-taxonomy-instrument.md`. *This is the
   `feedback_diagnostic_first` discipline — the lesson that cracked every hard bug
   in this lineage. We build S1 vs S2 first based on the measured split, not the
   comment's hint.*
1. **S1 — scale consensus only** (min-scale + Lipschitz, same-key edges). Expected:
   the bulk of bracket `skippedConflict` recovered (the dominant cause is scale
   divergence on shared rails). Gate: bracket conflict skips drop sharply,
   checkMesh-clean, no quality regression.
2. **S2 — canonical joint-path + assign/nSeg reconciliation.** Recovers the
   diverged-path conflicts. Gate: bracket conflict < 50; lbend family clean.
3. **S3 — convex-ridge seam** (ridge adjacency + shared seam). Gate: cow_ear
   `skippedConvexRidge` → 0, checkMesh-clean; if not fully clean, keep the skip
   fallback and hand the residual to I1.3.

Staging lets each lever be measured in isolation (diagnostic-first) and keeps every
intermediate build shippable.

## 6. Decisions for review (recommendations)

- **D1 — consensus seeding vs full rewrite. DECIDED 2026-06-16: SEEDING** (David:
  "the more robust, long-term option"). **Verified** the deciding fact: PASS 3 ring
  construction reads `RailData& rd = railsOf[railKey(col,c)]` (L798) and builds
  every ring from `rd.scale`/`rd.assign` — so `railsOf` is **already the single
  source of truth for construction**, not just validation. Seeding the negotiated
  record into `railsOf` therefore makes consensus authoritative for **both** passes
  with no change to the forensics-hardened validation/construction code. The "full
  rewrite" would add regression risk there for no architectural gain — rejected.
- **D2 — nSeg reconciliation on a shared rail.** *Recommend* **min-nSeg** at the
  seam (both members layer to the shallower common depth there) — robust, matches
  the "strictness over cleverness" lineage and the plan's "one negotiated (scale,
  assignment)". Alternative (depth-preserving, per-column nSeg with a shared
  prefix) is more complex; defer to I1.4 if coverage demands it.
- **D3 — convex-ridge scope I1.1 vs I1.3.** *Recommend* I1.1 builds the shared
  seam and **measures** on cow_ear: if checkMesh-clean → retire the skip now; else
  keep the skip as fallback and the fanned-corner geometry is I1.3. Avoids
  over-committing the hardest geometry to I1.1.
- **D4 — Lipschitz bound source.** *Recommend* the spec `expansionRatio` growth as
  the neighbour-scale bound (mirror `check_lipschitz`); calibrate on the bracket.

## 7. Risks & past mistakes to carry in

- **Verify the instrument first.** Measure recovery with the I4.2 comparator +
  `splitLayersReport` (now committed) and the calibrated checkMesh instrument —
  not by eye. Half a spike night was once a forensic-parser bug.
- **Determinism.** Union-find + min-scale + "most-cells" canonical path must be
  order-independent; tie-break by lowest column id / lowest point label so two runs
  are byte-identical (I4.3 fingerprint asserts it).
- **Grading interaction.** Lipschitz smoothing must not *raise* any scale (only
  lower toward the component min) or it re-introduces the height the column could
  not achieve — keep it a min-monotone smoother.
- **Never write a checkMesh-failing mesh.** The I4.1 self-check gate (committed)
  refuses a bad mesh regardless; consensus must turn skips into *clean* layers, and
  the matrix gates it.
- **GPL wall.** All of this is executor mechanics inside the GPL `splitLayers`; it
  stays behind the D7 process boundary. No planner intelligence moves here.

## 8. Fixtures (capture into the I4.2 matrix on the way)

- **bracket_corner** (captured bug): must FAIL on the pre-S1/S2 binary (conflict
  skips ~345, the rail-identity defect) and PASS after (conflict < 50,
  checkMesh-clean, nonOrtho ≤ 70).
- **cow_ear** (captured bug): must FAIL pre-S3 (`skippedConvexRidge` > 0 or skew)
  and PASS after (`skippedConvexRidge` → 0, checkMesh-clean, skew ≤ 4).
- lbend family + cavity: must stay green throughout (no regression).
