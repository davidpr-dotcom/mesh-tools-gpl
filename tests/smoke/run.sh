#!/usr/bin/env bash
# Smoke test: mesh the cavity case, run splitLayers in dry-run mode against
# the example spec, assert the structured diagnostics appear.
set -euo pipefail
cd "$(dirname "$0")/cavity"

blockMesh > log.blockMesh 2>&1

splitLayers -spec ../specs/example-layer-spec.json -dryRun > log.splitLayers 2>&1

grep -q "SPLITLAYERS|patch|movingWall" log.splitLayers
grep -q "SPLITLAYERS|result|dry-run-ok" log.splitLayers

echo "smoke: OK"
