#!/usr/bin/env bash
# I4.2 regression matrix — RUN half (GPL/WSL). Produces log.splitLayers +
# log.checkMesh for every fixture case; the JUDGE half is meshcore-side
# (meshcore.executors.regress.judge_matrix), so the GPL repo stays free of any
# proprietary dependency. WSL-only (OpenFOAM v2512 + splitLayers on PATH).
#
# Layout:  tests/regress/cases/<fixture-name>/run.sh   (one per manifest fixture)
#          each run.sh meshes its case, runs splitLayers (edge / optimizeSweeps 0),
#          and checkMesh, writing log.splitLayers + log.checkMesh into its dir.
#
# Usage:   bash tests/regress/run-all.sh
#          ZEFRA_REPO=/mnt/c/Users/traca/Desktop/fullCrossflow/Zefra \
#              bash tests/regress/run-all.sh        # also runs the meshcore judge
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
cases="$here/cases"
rc=0

[ -d "$cases" ] || { echo "no cases/ dir yet — capture fixtures first (see README.md)"; exit 2; }

for d in "$cases"/*/; do
    name="$(basename "$d")"
    if [ ! -x "$d/run.sh" ]; then
        echo "SKIP  $name (no run.sh)"; continue
    fi
    echo "RUN   $name"
    ( cd "$d" && bash run.sh ) > "$d/log.run" 2>&1 \
        || { echo "      tool error in $name (see log.run)"; rc=1; }
done

# JUDGE (optional, meshcore-side): parse the logs + compare against the manifest.
if [ -n "${ZEFRA_REPO:-}" ]; then
    echo "--- meshcore judge ---"
    PYTHONPATH="$ZEFRA_REPO" python3 - "$cases" <<'PY'
import sys
from meshcore.executors.regress import judge_matrix
m = judge_matrix(sys.argv[1])
for r in m["results"]:
    flag = "OK  " if r["ok"] else "FAIL"
    print(f"  {flag} {r['name']}: {r['metrics']}")
    for f in r["failures"]:
        print(f"        - {f}")
if m["missing"]:
    print("  MISSING:", ", ".join(m["missing"]))
print(f"matrix: {m['n']} judged, {m['n_failed']} failed, {m['n_missing']} missing")
sys.exit(0 if m["ok"] else 1)
PY
    rc=$(( rc | $? ))
fi

exit "$rc"
