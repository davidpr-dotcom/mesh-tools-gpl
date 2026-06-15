# splitLayers regression matrix (I4.2)

The **trust scaffold** that protects every coverage change (I1.x) from regressing.
Two halves:

- **RUN** (this dir, GPL/WSL): `run-all.sh` meshes each fixture, runs `splitLayers`
  (`-placement edge -optimizeSweeps 0`) and `checkMesh`, writing
  `log.splitLayers` + `log.checkMesh` per case.
- **JUDGE** (meshcore, proprietary): `meshcore.executors.regress.judge_matrix`
  parses those logs (`splitlayers_report` + `checkmesh`) and compares each against
  `meshcore/executors/contracts/regression_manifest.v1.json`. The GPL repo has no
  meshcore dependency; `run-all.sh` only invokes the judge when `ZEFRA_REPO` is set.

## Fixtures (manifest names → case dirs `cases/<name>/`)

Each case dir holds a `run.sh` that produces the two logs. Adapt
`tests/smoke/run.sh` (cavity, blockMesh) and `tests/smoke/run-lbend.sh`
(lbend, blockMesh + perturb) for:

- `lbend_clean`, `lbend_graded`, `lbend_opposing`, `lbend_perturbed`,
  `lbend_refineMesh` — the L-bend family (blockMesh variants).
- `cavity` — the existing smoke baseline.
- **`bracket_corner`** — CAPTURED BUG (2026-06-12 rail-identity). Capture from the
  WP 5.4 spike cfMesh base (`myCases/wp54_spike/.../concave_bracket`). Must FAIL on
  the pre-fix binary, PASS after I1.1.
- **`cow_ear`** — CAPTURED BUG (2026-06-13 `strips:opposed` convex ridge). Capture
  from the spike cow base; the I0.1 doc points at the freshly-built `split_dbg`.
  Must FAIL on the pre-fix binary, PASS after I1.1/I1.3 (`skippedConvexRidge → 0`).

## Capturing the two captured-bug fixtures (the must-fail check)

The DoD: each captured fixture **fails on the pre-fix binary and passes after**.

1. Build the **pre-fix** `splitLayers` (git checkout the commit before the
   I1.1/I1.3 fix), run the fixture → confirm the judge reports FAIL
   (bracket: conflict skips / nonOrtho; cow: `skippedConvexRidge > 0` or skew).
2. Build the **current** binary, run → confirm PASS.
3. Commit the case dir (mesh inputs + `run.sh`), not the logs (regenerated).

## Run

```bash
# RUN only (produce logs):
bash tests/regress/run-all.sh
# RUN + JUDGE:
ZEFRA_REPO=/mnt/c/Users/traca/Desktop/fullCrossflow/Zefra bash tests/regress/run-all.sh
```

Also fold the **AMR guard** (I0.2 `refineUnrefineCycle`, with the `hexMatcher`
filter so 6-face non-hexes don't crash `hexRef8`) into the matrix as a separate
check on the bracket + cow split meshes (zero destroyed layer cells).
