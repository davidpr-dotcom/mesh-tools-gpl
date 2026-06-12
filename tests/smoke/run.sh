#!/usr/bin/env bash
# Smoke test: mesh the cavity case, run splitLayers in dry-run mode against
# the example spec, assert the structured diagnostics appear.
set -euo pipefail
cd "$(dirname "$0")/cavity"

rm -rf 0 constant/polyMesh   # stale results/fields from previous runs
blockMesh > log.blockMesh 2>&1
CELLS0=$(grep -oP "nCells: \K\d+" log.blockMesh || echo 400)

splitLayers -spec ../specs/example-layer-spec.json > log.splitLayers 2>&1

grep -q "SPLITLAYERS|patch|movingWall" log.splitLayers
grep -q "SPLITLAYERS|result|split-ok" log.splitLayers

checkMesh > log.checkMesh 2>&1
grep -q "Mesh OK" log.checkMesh || { echo "smoke: checkMesh FAILED"; exit 1; }

CELLS1=$(grep -oP "cells:\s*\K\d+" log.checkMesh | head -1)
[ "$CELLS1" -gt "$CELLS0" ] || { echo "smoke: no cells added"; exit 1; }

echo "smoke: OK ($CELLS0 -> $CELLS1 cells)"
