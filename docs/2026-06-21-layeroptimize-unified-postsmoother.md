# layerOptimize — unified never-worse layer post-smoother (H2b design)

**Date:** 2026-06-21
**Status:** DESIGN — await approval; build is incremental, spike FIRST (prove-before-modify).
**Goal:** reduce the extrusion arm's internal nonOrtho (~75, passes checkMesh = polish) with a never-worse, BL-preserving post-insertion smoother — and make it the **one optimizer that serves both arms** (split already has it inline; this gives extrusion the same lever and unifies them).
**Builds on:** the committed split optimizer (localDescent + cornerScale + active-set + early-stop in `splitLayers.C`) — the algorithm is proven; this ports it to operate on an arbitrary polyMesh.

---

## 1. Why a standalone app

The split arm optimizes using its in-memory construction metadata (`inserted`/`layerOf`/`baseNewOf`/rails). The extrusion arm's output is just a **polyMesh + a `layerCells` cellSet** (snappy's added cells, already emitted by `snappy_layers.compute_layer_cells` → `layerZoneExtract`) with **no such metadata**. So a standalone GPL app — `layerOptimize` — that reads `(polyMesh, layerCells cellSet, wall patch)` and reconstructs what the descent needs is the clean way to give extrusion the proven optimizer. It then also works as a post-step on *any* layer mesh (split included) → the unified post-smoother.

## 2. The one new piece — prism-column detection from topology (SPIKE THIS FIRST)

Everything else is ported + proven; the only new algorithm is recovering the column structure the descent + y⁺-pinning need:

- For each **wall-patch face** `wf` (the base of a column), its owner cell `c0` is the first layer cell.
- **Walk outward:** from layer cell `c_k`, take the face *opposite* the incoming face (for a hex prism, the face sharing no vertices with the incoming face); its neighbour is `c_{k+1}`. Continue while the neighbour is in `layerCells`; stop at the first non-layer cell (the bulk) or a boundary.
- This yields, per column: cells `[c0..c_{n-1}]`, ring faces `[f0=wf, f1, …, fn]` where `fk` points = the **layer-k points**. `f0` (wall) points are FIXED; `f1..fn` points are the movable layer points, each tagged with its layer index `k`, its column, and the **wall point directly below it** (the `f0` corner connected by the prism's vertical edge) = its y⁺ base.
- **Non-hex columns** (snappy near features) are skipped (same posture as splitLayers' non-hex skip). v1 targets the hex prism stacks that carry the bulk of the layer faces.

The spike proves this recovers correct columns (count ≈ wall faces; per-column `k` monotone; wall-point mapping correct) on a real snappy internal mesh **before** we wire the descent to it.

## 3. Optimize (ported, proven)

Once each movable layer point has `{wall point (y⁺ base), k, column}`, run the **exact** committed descent:
- **localDescent** — active-set faceQ never-worse descent (T = max(30, 0.7·curMax), early-stop tol 1e-4), y⁺-pinned by renormalizing each move to keep the wall-offset magnitude `|p − wallPoint|`.
- **cornerScale** — per-column radial re-spacing (k=0 pinned), compress + mild expand, never-worse.
- **stall-driven pipeline** — both to stall, repeat until no round-progress.

`faceQ`/`localJ`/`maxLayerNonOrtho`/the descent/active-set/cornerScale are lifted from `splitLayers.C`. v1: **copy** them into `layerOptimize` to avoid touching the committed splitLayers; once green, **extract to a shared header** (`layerOptCore.H`) so both apps share one implementation (no drift).

## 4. Invariants (unchanged from the split optimizer)

- **Never-worse** — global max nonOrtho + skew monotone non-increasing (the colQ/localJ gate scans pointCells of every moved point).
- **y⁺ preserved** — wall-offset magnitude held per layer point (k=0 / wall points never move); this is the BL-preservation that unconstrained centroidal smoothers (smoothMesh's trap) lack.
- **No tangle** — never-worse rejects inverting moves; cornerScale compresses by default.
- **Deterministic** — same input → same output.

## 5. Incremental plan (each step gated before the next)

1. **Spike the column detection** (§2): a small `-detectOnly` mode (or a probe) that walks the stacks on a snappy internal mesh (e.g. straight_pipe/u_bend at ~75) and emits `LAYEROPT|columns|N|layerPts|M|nonHexSkipped|K` + a sanity dump. Confirm columns ≈ wall faces, monotone k, correct wall mapping.
2. **Build `layerOptimize` v1** (localDescent only) on the detected points; run on the spike mesh.
3. **Validate**: internal nonOrtho drops, **first-layer thickness (y⁺) preserved**, never-worse, coverage unchanged (it moves points, not layers).
4. **Add cornerScale** if v1 under-shoots.
5. **Wire** into `snappy_layers.add_layers` as an optional post-step (run `layerOptimize` after addLayers, re-checkMesh) + re-run the gate. Apply to the split arm too if it beats the inline path (the unification).

## 6. Scope / honesty

Internal extrusion nonOrtho ~75 **passes checkMesh today** — this is solver-conditioning polish, not a gate fix. The real prize is the **unified, proven, never-worse post-smoother** as durable infrastructure. v1 is hex-prism columns on internal extrusion; external already won via H1 params. If the detection proves fragile on snappy's non-hex feature cells, we fall back to evaluating smoothMesh (H2a) — but the spike tells us that cheaply, first.

## 7. Validation targets
straight_pipe_00, u_bend_00, annular_duct_00 (internal extrusion, ~75). Metric: maxNonOrtho ↓, first-layer thickness within ~1% of pre-smooth (y⁺), checkMesh-clean, never-worse vs the un-smoothed layered mesh.
