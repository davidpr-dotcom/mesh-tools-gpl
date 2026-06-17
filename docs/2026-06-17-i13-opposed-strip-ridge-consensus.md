# splitLayers: I1.3 — opposed-strip convex-ridge consensus (2026-06-17)

**Track:** Phase 5 / boundary-layer engine. The real fix for the I0.1 cow-ear
defect, now the dominant blocker once the geometric convex skip is retired.
**Status:** design, pre-implementation. Needs David's approval before C++.
**Replaces:** the `-convexRidgeAngle` geometric skip (a mask) + the orchestrator
demote (a fallback that loses the layers).

## The defect (root-caused, WSL 2026-06-17)

Retiring the convex skip exposed the manifold's bulk failure: **17790 open cells,
all `nCol|8` (full, NOT gap-capped), all tagged `strips:seg0:cuts8:ocNNN:oseg0`**
— two columns each meeting strip-to-strip at a convex ridge. The cow's 58 open
cells are the same class. It is independent of the gap-cap (those columns reach
full depth, so they're never seen as gaps).

## Why opposed strips don't close (grounded in the code)

- Rings are created per `railKey = edge(rails[c][0], rails[c][1])` = (wall point,
  **first interior joint**) and stored in `railsOf`. Two columns that produce the
  *same* railKey **share** the same ring labels.
- `withForeignRings` only shares rings along **rail segment edges** (`segRings`),
  i.e. along a rail, not along a wall edge.
- At a **flat** shared wall point P, adjacent columns A,B grow rails from P into the
  *same* cell → *same* first interior → *same* railKey → **shared rings** → the
  side face on the shared wall edge is single-valued → closes. (This is why flat
  regions are fine.)
- At a **convex ridge**, A's rail from P bends along A's wall and B's along B's →
  *different* first interiors → *different* railKeys → **separate, divergent ring
  sets**. The shared side face (built once, with one column's rings) then leaves
  the other column's layer cells unclosed → the `strips:opposed` open/skew faces.

So the defect is **per-wall-point**: at a convex-ridge point, the columns around it
each invent their own divergent ring column instead of sharing one.

## The fix — consensus rings at ridge points

At a convex-ridge wall point P, create **one** ring column, placed along the
**smoothed inward normal** at P (`patchNormals` already computes this — it is the
bisector of the faces around P), and have **every** column touching P reuse those
ring labels. The strips around P then meet on a single shared ring column → close.

Concretely:

1. **Ridge-point set.** `smoothedInwardNormals(pp, …, featAngle, nFeat)` already
   detects feature points (the `featurePoints` count). Expose/derive the set of
   convex feature points (sharp + convex, distinguished from concave by the
   face-centre test already used in the old convex-skip block, L597-628).
2. **Consensus ring registry.** `Map<label, labelList> ridgeRings` (wall point P
   → ordered ring labels along the smoothed normal: `P + n̂_P · cumHeights[k]`,
   `k=0..nCol-1`). Create lazily the first time a column needs P's rings; all
   later columns at P reuse them. (Mirrors the per-segment `segRings` reuse that
   S2b leans on, but keyed by **wall point** for ridges.)
3. **Ring creation (PASS 3).** When placing a rail's ring `k` whose **base wall
   point** is a ridge point, take it from `ridgeRings[P]` instead of creating a
   per-railKey point. Non-ridge (flat/interior) points keep the existing rail
   placement → flat regions byte-identical.
4. **Strip pass (section 5).** Unchanged — it already consumes whatever ring
   labels the cell carries; with `ridgeRings` shared, the opposed strips now use
   identical points and close. `withForeignRings` likewise picks up the shared
   labels on the wall edge if we register the wall-edge ring set alongside.

This makes the convex ridge a **single consensus column fan** (the standard
commercial approach) rather than N divergent per-column rails.

## Staging (each WSL-built + verified; diagnostics stay on)

- **I3.1 — ridge-point detection, emit-only.** Build the convex-ridge point set,
  emit `SPLITLAYERS|ridgePts|<n>`. Verify it matches the open-cell clusters
  (manifold high, ahmed 0). Proves the detector before any construction change.
- **I3.2 — consensus ring registry + ring creation.** `ridgeRings`, lazily filled
  along the smoothed normal; PASS-3 ring placement reads it at ridge base points.
  Verify: manifold/cow opposed-strip open cells → ~0; **ahmed/bracket byte-
  identical** (no convex ridges → no ridge points → unchanged) — the mechanical
  no-external-regression gate.
- **I3.3 — wall-edge ring sharing for `withForeignRings`** if any residual edges
  remain; tune.

## Verification / DoD

- manifold + cow: `strips:opposed` open cells → 0; full layers present at convex
  ridges (no demote, no skip); checkMesh clean (or orchestrator-clean for genuine
  residual skew only).
- ahmed/bracket: **byte-identical** (no convex ridges).
- With this landed, the `-convexRidgeAngle` skip can stay retired (default 179)
  and all four geometries are clean with full layers — the no-trade-off goal.
- The OPENCELL / badFace / unhandled / OOB diagnostics remain as guards.

## Risks

- **Touches ring creation** (PASS 3) — the core. Mitigated: ridge points are a
  small set (feature points only), flat/interior placement is unchanged
  (byte-identical externals), and the open-cell diagnostic catches any new
  non-closure at the source.
- **Ridge vs concave discrimination** must be right (concave corners must NOT get
  the convex consensus); reuse the existing face-centre convex test (L614).
- **Determinism**: `ridgeRings` filled in column order, keyed by point → order-
  independent; assert via the report fingerprint.
- **Overlap with I1.1 S2b** (corner bisectors): related but distinct — S2b is the
  reentrant/divergent-suffix prefix case; this is the convex-ridge fan. Keep both;
  they compose (both make adjacent columns share ring labels).
