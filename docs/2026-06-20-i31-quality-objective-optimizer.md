# I3.1 — quality-objective layer optimizer (design)

**Date:** 2026-06-20
**Status:** DESIGN — await approval before C++ (David builds/runs on WSL).
**File:** `applications/splitLayers/splitLayers.C`, the offset-vector smoothing pass (~L2177–2269).
**Plan ref:** `Zefra/meshcore/docs/2026-06-15-phase5-layer-engine-plan.md` §I3.

---

## 1. Why (measured, not assumed)

The Stage-C nonOrtho attribution (`Zefra/tools/2026-06-20-nonortho-attribution.py`, base-vs-layered checkMesh) localized the standing maxNonOrtho:

- **Smooth bodies / extrusion** (prolate, proper mesh): base 18.7 → layered 26.6 — clean. The raw-corpus ~74 was a single-patch mesh artifact. **No target here.**
- **Sharp-feature bodies / split** (ahmed): base **62.0** → layered **80.0** (Δ +18), skewness **1.37 → 3.94**. The split layers genuinely *introduce* nonOrtho and skew. **This is I3.1's target.**

So I3.1 applies to split-routed sharp-feature geometry (the MVP external-aero bodies — cars, wings). It can reduce 80 toward the **base floor ~62** (the cartesian cut at sharp edges, a separate base-mesh concern), never below it. The 80 is checkMesh-clean today, so this is quality polish (better solver conditioning), not a gate fix.

## 2. What's there now (the blind smoother)

`if (nOptSweeps > 0)` builds, per rail, the inserted layer points indexed by layer `k` (`inserted`, `layerOf`, `baseNewOf`), then for each sweep:

```
v        = 0.5*vOwn + 0.5*(mean same-layer neighbour offset)   // vOwn = newPts[np]-base
target   = base + v renormalized to |vOwn|                     // out-of-plane (geomD<0) frozen
newPts[np] = target                                            // moved UNCONDITIONALLY
```

It is a tangential Laplacian on the offset vectors with **no quality test** — on warped/sharp walls the averaged direction can tilt a face worse, so full sweeps regressed bracket 43.9 → 78.9. Hence `sweeps=0` ships.

## 3. The rebuild — accept-only-if-better

Keep the smoothing *candidate* (the existing `target`); gate the move on a measured local objective.

**Local objective.** For an inserted point `np`, over its incident INTERNAL faces `f` (owner `o`, neighbour `n`):

```
J(np) = max_f  max( nonOrtho(f) / NONORTHO_REF , skewness(f) / SKEW_REF )
        nonOrtho(f) = angle( C[n] − C[o] , Sf )           // checkMesh's definition
        NONORTHO_REF = 70°,  SKEW_REF = 4.0                // checkMesh severe thresholds
```

Normalizing by the checkMesh limits lets one objective guard **both** nonOrtho and skew (ahmed's skew also jumped). Boundary faces are skipped.

**Accept rule (Gauss-Seidel, in place).** Per sweep, per inserted `np`:
1. `Jcur  = J(np)` at the current position.
2. `Jcand = J(np)` with `np` provisionally at `target` (local geometry recomputed — see §4).
3. Move iff `Jcand <= Jcur + EPS`; else leave `np` put. Out-of-plane components stay frozen (`geomD<0`).

**Freeze near skips.** Skip any `np` whose incident cells include a skipped/forced-skip column (those interfaces are fragile — the same reason the recovery froze them). Carry the skip-column set already tracked in PASS-2/quality-retry.

**Why this is never-worse (the DoD guarantee).** Before a move, every face `np` touches is ≤ `Jcur ≤ global max`. We accept only if after the move all of `np`'s faces are ≤ `Jcur` (that's what `Jcand ≤ Jcur` means, since `J` is the max over those faces). Faces not incident to `np` are unchanged. So no face ever exceeds the pre-move global max ⇒ **global max nonOrtho is monotonically non-increasing across the whole pass.** It can only improve or hold — safe to default-on.

## 4. Implementation core — the local recompute

The gate must use checkMesh's metric, so `Jcand` recomputes the **exact** geometry of the bounded set affected by moving `np`:
- affected cells = `mesh.pointCells()[np]`; affected faces = `mesh.pointFaces()[np]` (∪ the faces of those cells needed for cell centroids).
- recompute those faces' area-normals/centres and those cells' volume-centroids with `np` at `target` (OpenFOAM `primitiveMeshTools`/`face::areaNormal`+`face::centre`, cell volume-centroid from its faces), then `nonOrtho`/`skewness` per incident internal face.
- everything else in the mesh is untouched, so this is O(incident faces) per point.

This local exact-recompute is the careful part and will WSL-iterate (getting the cell volume-centroid update right for a single moved point). Emit `SPLITLAYERS|optimize|sweeps|N|points|M|moved|K|maxNonOrthoBefore|x|maxNonOrthoAfter|y` so the effect is measured each run.

## 5. DoD

- **ahmed (split): max nonOrtho strictly DECREASES** (80 → toward the 62 base floor) with `sweeps>0`; skew does not increase. checkMesh-clean.
- **never-worse on every geom:** bracket / cow / wing / manifold max nonOrtho with `sweeps>0` is ≤ the `sweeps=0` value (the old regression — bracket 43.9→78.9 — must NOT recur).
- **deterministic:** two runs byte-identical (I4.3 fingerprint).
- Once the DoD holds across the spike suite, **flip the default** `optimizeSweeps` 0 → (e.g.) 10; until then it stays 0 and is opt-in.

## 6. Scope / caveats

- Helps the split-routed sharp-feature subset only (where the layers add nonOrtho). Smooth/extrusion geoms have no target (prolate clean at 26.6).
- Cannot go below the **base edge-cut floor** (ahmed ~62) — that's a separate, lower-priority base-mesh lever (improveMeshQuality / feature-edge cut handling), and it's checkMesh-clean.
- **I3.2 (follow-on):** retire `-placement normal` for the planner's smoothed direction field; orthogonal to this objective-gated optimizer and deferred.

## 7. Validation

Author the pass → David `wmake` → `spike-split-eval` (ahmed/bracket/cow + a wing) at `--optimizeSweeps 10 --checkmesh` vs `--optimizeSweeps 0`; judge by the emitted before/after maxNonOrtho per geom + checkMesh-clean + the never-worse comparison. Gate on the `bracket_corner`/`cow_ear` I4.2 fixtures (no regression).
