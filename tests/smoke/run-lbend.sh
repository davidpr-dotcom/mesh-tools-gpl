#!/usr/bin/env bash
# L-bend multi-row reproducer (diagnostic, not pass/fail):
# uniform 5 mm cells, spec forces 3-cell chains, reentrant corner included.
set -euo pipefail
cd "$(dirname "$0")/lbend"

rm -rf 0 constant/polyMesh
blockMesh > log.blockMesh 2>&1
python3 perturb.py
splitLayers -spec layer-spec.json -placement edge -optimizeSweeps 0 \
    > log.splitLayers 2>&1 || true
grep "SPLITLAYERS|" log.splitLayers
checkMesh > log.checkMesh 2>&1 || true
sed -n '/Checking geometry/,$p' log.checkMesh
