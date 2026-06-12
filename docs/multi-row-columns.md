# splitLayers: multi-row columns (design, 2026-06-12)

## Problem

Split-based layers carve within the wall-adjacent cell, so achievable layer
heights are bounded by wall-cell height. cfMesh trims wall cells (bracket:
~1–3 mm vs a 5.6 mm stack) → 65% first-cell error. Extrusion has no such
bound. Fix: extend the split column through a CHAIN of stacked cells until
the stack (plus a healthy remainder) fits.

## Model

A column is a chain of hex cells walked from the wall face through
successive opposing faces. Its four corner point-paths are RAILS
(polylines, one joint per cell interface). Ring points are placed along
rails at arc-length positions `scale * cumHeight_k`, where
`scale = min(1, railLen / (total + lastLayer*er))` — the remainder always
continues the geometric progression.

Topology per column (m cells, n layers): all m cells removed, n+1 stack
cells added; the m−1 interface faces are REMOVED (absorbed); wall face →
stack[0]; top face → remainder; each original side face is partitioned by
the rings assigned to its segment (point substitution preserves winding);
rail joints remain as polygon vertices for lateral neighbours (no hanging
points); lateral faces of unsplit cells get ring points inserted per
segment edge.

## Consistency rules (v1 — strictness over cleverness)

Chains stop at: boundary, non-hex/non-quad-faced cell, any wall-adjacent
cell of a target patch, any cell already claimed by another chain, rail
continuation failure.

Skip classes (reported per patch; skip rate is a primary spike metric):
- `skippedCorner` — cell owns >1 target patch face (as before)
- `skippedNonHex` — wall cell not splittable (as before)
- `skippedConflict` — a shared rail already subdivided with a different
  scale/segment-count (neighbouring chains of unequal length); rings are
  a property of the rail, first column wins
- `skippedWarp` — the four rails of a column disagree on which segment a
  ring falls in (heavily warped chains)
- `skippedCross` — a side face is claimed by another column with a
  different rail pair (perpendicular chains crossing, e.g. internal-flow
  cases near corners)

All conflict detection happens BEFORE point creation (identification →
conflict pre-passes → ring creation → topology), so a skipped column
leaves zero residue.

Placement: rings follow the rail polyline (arc-length). `-placement
normal` applies only to single-cell columns where it is proven; multi-row
normal placement deferred. Optimizer unchanged (offset-vector smoothing).

## Knowingly deferred

- Chains through octree level transitions (2:1 interfaces): chains stop.
- Relaxing skippedConflict via per-region scale consensus → planner-side
  thickness field (WP 5.1 item 3) is the real solution.
