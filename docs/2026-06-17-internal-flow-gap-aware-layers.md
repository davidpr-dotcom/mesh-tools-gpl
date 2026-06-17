# splitLayers: gap-aware layers for internal flows (I-INT, 2026-06-17)

**Track:** Phase 5 / boundary-layer engine — internal-flow quality.
**Status:** design, pre-implementation. Needs David's approval before C++.
**Reads with:** the I1.1 corner docs; this is the *internal* counterpart.

## Problem (grounded in the WSL attribution, 2026-06-16)

With the geometric convex skip retired (`convexRidgeAngle 179`, orchestrator owns
quality), the three **external** geometries are clean (ahmed/bracket/cow; bracket
recovery 74→131). The **internal** geometry (`manifold_gap`, a thin gap) fails
checkMesh hard: open cells 29512, nonOrtho 35172, wrongOriented 18551, skew 15740.

`tools/2026-06-16-attribute-manifold-defect.py` maps the offenders to the
constructs that built them:

- `strips` 33713, `wall+rings` 14756 (the two dominant families)
- `wall+rings` 14756 = **wall faces wrong-oriented** → the first layer cell is on
  the wrong side → **layers inverting**.
- `strips:seg3/seg10:cuts0:plain` 1218 each = **long-multiRow** side faces.
- the opposed-ridge (`oc…:oseg…`) and top-meeting tags are 2 faces each — NOT the
  cause.

Combined with `minScale 0.0306 / meanScale 0.179`: the spec's **8 layers**
(`y1 9.757e-05`, er 1.2 → total ≈ 1.6 mm) are being crammed into **sub-millimetre
gaps**. Layers from the two opposing walls grow toward each other and **overlap /
invert** through the gap centre; cells go degenerate. The quality-retry
orchestrator structurally **cannot** fix this — demoting one wall's columns just
opens the other's, so it cascades and exhausts (26839 demoted, still failing).

**Conclusion:** this is not a per-column quality defect (orchestrator's domain) —
it is a *layer-budget* defect. A gap that is too thin physically should not carry
8 layers. The fix belongs at **layer sizing**, and it is the case where internal
flow genuinely differs from external.

## The distinguishing principle (internal vs external)

For a wall face, ray-cast the **inward** normal and measure the distance `d` to the
nearest **opposing** wall (another face of the same layered patch whose normal
opposes ours). Then:

- **External** (open domain): no opposing wall within reach → `d = ∞` → no
  constraint → full layer stack (today's behaviour, unchanged).
- **Internal** (gap): `d` is the local gap width → the layers from both walls must
  share `d`; each wall gets at most ~`d/2` (minus a core margin).

This single measurement *is* the internal/external discriminator — no flow-type
flag needed; the geometry tells us. (External faces simply never find a near
opposing wall.)

## Design

### 1. Gap detection (`d` per wall face)

Build an `indexedOctree<treeDataFace>` (or `triSurfaceSearch`) over the layered
patch's faces. For each wall face `f`: cast a ray from its centre `C_f` along the
**inward** smoothed normal `-n_f` (we already compute `patchNormals` at L852); the
first patch face hit at distance `d_f` whose normal **opposes** `n_f`
(`n_f · n_hit < -cosOpp`, e.g. cosOpp = cos(60°)) is the opposing wall. No hit, or
the nearest hit is far (`d_f > D_external`, a multiple of the full stack height) →
`d_f = ∞` (external). Cost O(N log N); the octree is built once.

Emit `SPLITLAYERS|gap|<patch>|internalFaces|<count>|minGap|<x>` (diagnostic).

### 2. Per-column layer-count cap (the budget)

Full stack height for `m` layers: `T(m) = y1·(er^m − 1)/(er − 1)`. Choose, per
column, the largest `n_col ≤ nLayers` with `T(n_col) ≤ d_f/2 − margin`
(`margin` = a fraction of `d_f`, leaves a clean core band the two stacks meet in).
`d_f = ∞` ⇒ `n_col = nLayers` (external, unchanged). Very thin gap ⇒ `n_col` small,
possibly 0 (no layers — a clean, deterministic skip, far better than the inverted
cells the convex crutch was hiding).

This is the **only new sizing input**; `col.scale = min(1, minRailLen/needed)`
still applies, now over `n_col` layers.

### 3. Variable-`n` topology (the real work)

Today every column builds `nLayers+1` stack cells and the ring/strip/topology
assume a **uniform** `n` across neighbours. With per-column `n_col` we need
neighbouring columns of **different** counts to mesh — standard "layer
termination" (snappyHexMesh does this):

- **Stack:** `n_col + 1` cells per column (n_col layers + remainder). (L1197-1209.)
- **Ring faces (L1328-1338):** loop `k = 0 .. n_col-1` — unchanged but bounded by
  the column's own `n_col`.
- **Side/strip faces (L1459+, the crux):** for a face between column A (n_A) and
  neighbour B (n_B), sub-segment `s`: A's side maps to `stackA[ringsBeforeA+s]`,
  B's side to `stackB[ringsBeforeB+s]`. When `s` exceeds B's layer count, B's layer
  doesn't exist → map to **B's remainder** `stackB[n_B]` (A's upper layers
  *terminate* against B's core). This is exactly the clamp the OOB instrument
  already does (`chk` → size-1 = remainder); the change is to make it the
  **intended** termination (clamp to `min(idx, n_B)` → remainder), not a
  bug-guard. Geometrically the terminating face is a thin-to-core transition
  (some non-ortho, acceptable + the orchestrator can still demote true outliers).
- **Top face:** unchanged (the defer fix stands).

### 4. Where it lives

PASS 1 (where `col.scale` and the per-column data are computed) gains `col.nCol`
(the capped count) from the gap field. The spec's `nLayers` becomes the **cap**,
not a fixed value. The octree gap pass runs once before PASS 1's per-column loop.

## Staging (each step WSL-built + verified, behind the existing diagnostics)

- **G1 — gap field, emit-only.** Build the octree, compute `d_f`, emit the `gap`
  diagnostic. No behaviour change. Verify: external geoms report `internalFaces 0`
  (or few); manifold reports a large internal count + a sub-mm `minGap`. This
  *proves the discriminator* before we act on it.
- **G2 — count cap.** Apply `n_col` from `d_f`. Variable-`n` stacks + ring faces.
  Strip-pass termination (step 3). Verify: manifold checkMesh ≪ (open cells → ~0),
  external geoms **byte-identical** (their `n_col == nLayers`).
- **G3 — tune the margin / core band** so the manifold is clean (or
  orchestrator-clean) with the *most* layers that fit; confirm y⁺ is sensible.

## Verification / DoD

- Manifold: checkMesh clean (or orchestrator converges, not exhausts); layers
  present where the gap allows, gracefully fewer/none where it can't.
- ahmed/bracket/cow: **no regression** (byte-identical — they're external,
  `n_col == nLayers`). This is the "no trade-off" guarantee, checked mechanically.
- The I4.2 regression matrix + the quality-retry drivers stay green.
- Deterministic (the octree ray-cast + the budget are order-independent).

## Risks

- **Variable-`n` termination** is the riskiest change (touches the strip/topology
  core). Mitigations: the existing `badFace`/`unhandled`/OOB diagnostics stay on
  and will catch any dangling/over-run face at the source; G1 proves the
  discriminator before any topology change; external geoms are a byte-identical
  regression gate.
- **Ray-cast robustness** (grazing/curved gaps): use the smoothed normal + an
  opposing-normal test (`< -cosOpp`) so we only treat genuine facing walls as a
  gap; ambiguous → treat as external (full layers, today's behaviour) = safe.
- **Margin tuning** (G3) is a knob, but bounded and diagnostic-driven, and the
  orchestrator backstops residual quality — unlike the convex angle, a wrong
  margin degrades gracefully (fewer layers), never inverts.

## Out of scope (separate tracks)

The warped corner-bisector columns (I1.4, task 18) and the bracket's remaining
~91 warp aborts are unrelated to the gap budget and stay on their own track.
