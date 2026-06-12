#!/usr/bin/env python3
"""Deterministically warp the lbend mesh points (~20% of local cell size)
to emulate cfMesh's snapped/trimmed boundary cells. Keeps the outer hull:
points on the domain boundary planes only slide WITHIN their plane."""

import math
import re
from pathlib import Path

PM = Path(__file__).parent / "constant" / "polyMesh"
txt = (PM / "points").read_text()

dict_end = txt.find("}") + 1            # end of FoamFile header dict
body = txt[txt.find("(", dict_end) + 1:]

pts = [[float(x) for x in m.group(1).split()]
       for m in re.finditer(r"\(([-0-9.eE+ ]+)\)", body)]
header = txt[:dict_end]

A = 0.0009   # amplitude [m] ~ 20% of the finer cells


def warped(p):
    x, y, z = p
    dx = A * math.sin(417.0 * y + 93.0 * z + 1.3)
    dy = A * math.sin(389.0 * x + 71.0 * z + 2.1)
    dz = A * math.sin(403.0 * x + 59.0 * y + 0.7)
    # boundary planes: zero the normal component of the displacement
    eps = 1e-9
    if abs(x) < eps or abs(x - 0.06) < eps or abs(x - 0.02) < eps and y > 0.02 - eps:
        dx = 0.0
    if abs(y) < eps or abs(y - 0.06) < eps or abs(y - 0.02) < eps and x > 0.02 - eps:
        dy = 0.0
    if abs(z) < eps or abs(z - 0.02) < eps:
        dz = 0.0
    return (x + dx, y + dy, z + dz)


out = [header + "\n\n"]
out.append(f"{len(pts)}\n(\n")
for p in pts:
    w = warped(p)
    out.append(f"({w[0]} {w[1]} {w[2]})\n")
out.append(")\n")
(PM / "points").write_text("".join(out))
print(f"perturbed {len(pts)} points (A={A*1e3:.2f} mm)")
