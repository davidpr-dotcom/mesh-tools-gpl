# Optimize-phase dispatch — multiple optimizers, engine-selected (design)

**Date:** 2026-06-20
**Status:** DESIGN — await approval before code (David builds/runs on WSL).
**Idea (David):** make several layer optimizers, each suited to a different defect, merge them into one optimize phase, and let the mesh engine pick the one(s) that best fit each case.
**Builds on:** I3.1 v2 (committed) — the quality-objective local descent, now the first registered optimizer; and the WORSTFACE locator (committed) — the seed of the dispatch census.

---

## 1. Why this is the right shape

It is the **same pattern that already works for the layer engine**: split-vs-extrusion is "two engines, a coverage-aware dispatcher picks the one that suits the geometry." The optimize phase has the same structure — different defect classes want different algorithms — so the same architecture applies:

- a **registry** of optimizers, each with an *applicability* test and a *never-worse* guarantee;
- a **census** of the built layer mesh (what's limiting quality, and where);
- a **dispatcher** that selects + composes the best-suited optimizer(s) for this case.

The universal invariant that makes composition safe: **every optimizer is never-worse** (max nonOrtho AND max skew monotonically non-increasing — the I3.1 v2 property). So the engine can run several in sequence and quality only ever holds or improves.

## 2. Architecture

**Optimizer interface** (each registered optimizer):
- `name`.
- `targets` — the defect class (e.g. distributed-nonOrtho, corner wall+rings, poor normals, skew).
- `applies(census) -> score` — from the census, how much of the limiting defect this optimizer addresses (0 = not its job).
- `apply(mesh, region)` — run it, **gated never-worse** (the shared contract); emit before/after.

**Census** (the dispatch input — generalize the WORSTFACE locator): classify the worst-quality faces by *construct + cause + distribution*:
- construct tag (ridge/strip vs `wall+rings` vs lateral) via the faceMap+faceTagOf path already added;
- cause bucket (sharp-ridge / tight-corner-wall+rings / poor-normal / skew-dominated / non-hex);
- distribution (localized cluster of columns vs broadly distributed) — ahmed was a 4-column cluster, bracket is distributed.

Emit it as `SPLITLAYERS|census|...` so the dispatch decision is measured + logged (mirrors the dispatch dry-run line).

**Dispatcher** (in-tool — "the engine decides"): score each registered optimizer's `applies(census)`, run the applicable ones in priority order, re-censusing between (cheap), each never-worse. One merged phase; no planner round-trip.

## 3. Initial registry

| optimizer | targets | suits | status |
|-----------|---------|-------|--------|
| `localDescent` | distributed moderate nonOrtho/skew | bracket (45→41) | **DONE (I3.1 v2)** |
| `cornerScale` | `wall+rings` nonOrtho at tight corners | ahmed rear 3-way corners (~80) | NEW — validate |
| `directionField` (I3.2) | poor-normal-driven nonOrtho | corner/feature normals | LATER |
| *(future)* `globalSmooth`, `skewTargeted` | coupled / skew-led defects | — | OPEN |

`cornerScale` is the first real test of the dispatch: at the localized high-nonOrtho `wall+rings` columns (the census flags them), reduce the layer SCALE there (thinner rings ⇒ ring centroids hug the wall normal ⇒ lower nonOrtho), gated never-worse, accepting a local y⁺ trade-off only where it strictly helps. If it moves ahmed's ~80 where `localDescent` could not, the dispatch has earned its keep.

## 4. Invariants + DoD

- **Never-worse, composed:** each optimizer is gated; the phase is monotone on max nonOrtho + max skew. The old blind smoother's regression (bracket 43.9→78.9) can never recur.
- **Deterministic:** byte-identical across two runs (I4.3 fingerprint).
- **DoD:** the dispatched phase is **≥ the best single optimizer on every spike geom**, never-worse everywhere, deterministic, and `cornerScale` strictly improves at least one geom (ahmed) that `localDescent` left flat. Then default-on.

## 5. Incremental plan

1. **Refactor** the committed optimize pass into `localDescent` behind the optimizer interface; the dispatcher with one optimizer == today's behavior (byte-identical regression check).
2. **Census**: promote WORSTFACE into the structured census (construct + cause + cluster), emit `SPLITLAYERS|census|...`.
3. **`cornerScale`** + its `applies` (localized `wall+rings` clusters) → validate it reduces ahmed's ~80; confirm never-worse + no regression elsewhere.
4. **I3.2 `directionField`** only if a geom's census points to poor normals.

## 6. Honest scope

This session proved the standing quality gaps are largely **clean or geometrically irreducible** (prolate 26.6 clean; ahmed ~80 checkMesh-clean; bracket already passes). So the dispatch's *immediate* fix is modest — its real value is **architectural**: extensibility (new optimizers slot in behind one interface), per-case best-fit, and the same proven dispatch pattern as the layer engine. Recommend building it as the durable Stage-C+ optimize architecture, with `cornerScale` as the concrete proof that dispatch beats any single optimizer — not as a fix for a pressing defect (there isn't one today).

## 7. Implementation status (2026-06-20)

- **localDescent** + **active-set gate** (descend only points incident to a ≥T face, `T = max(30, 0.7·curMax)`, rebuilt each sweep) + **early-stop** (relative tol `1e-4`). Active-set proven equivalent on bracket (38.4586 exact, 6× fewer moves); cow `--sweeps 10` went 93s→33s — the descent was the bottleneck, so the build-once cache was dropped.
- **Survey corrected the targets:** cow = `wall+rings` (80.84), ahmed = `strips`/ridge (83). So cornerScale's real target is **cow**, not ahmed.
- **cornerScale (NEW):** per hot column (max face ≥T), re-space the OUTER layers (k≥1) by one compress factor `g ∈ {0.6,0.75,0.9}` while PINNING k=0 — so y⁺ is preserved *by construction* (no y⁺ floor needed). Gated never-worse via `colQ` over the column points' pointCells faces ⇒ global max non-increasing, same guarantee as localDescent. v1 is compress-only so it cannot tangle.
- **Dispatcher = stall-driven pipeline** (the "combine when one stalls" idea): each round runs localDescent→stall then cornerScale→stall; stop when a full round makes no global progress (`OPT_MAX_ROUNDS=4`). The census/`worstTag` is now **diagnostic** (the pipeline runs every lever and the census tells us which helped) — so it needs no tags and works in production without `--debug`.
- **Validation = cow.** Open question, tested not assumed: does the pipeline beat localDescent-alone's 80.84? cornerScale is a *hypothesis* (radial scale may or may not relieve `wall+rings`); never-worse makes the test safe. If it stalls too, scale isn't cow's lever → `directionField` (I3.2) is next.
