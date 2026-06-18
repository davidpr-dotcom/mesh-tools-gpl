# mesh-tools-gpl

OpenFOAM/cfMesh-linked command-line utilities supporting the Zefra meshing
engine. **GPL-3.0** — every binary in this repository compiles against
OpenFOAM (and, where noted, cfMesh) libraries and is published under their
license terms from day one (engine plan decision D10).

## What this repo is — and is not

These tools contain **zero proprietary intelligence**. They are mechanical
executors: each reads a versioned, schema-validated spec file (see
`etc/contracts/`), performs a mesh operation with OpenFOAM library calls, and
writes results plus structured diagnostics. All planning logic (sizing,
layer specs, retry policy) lives in a separate proprietary engine that
communicates with these tools **only via files across a process boundary**
(decision D7). No linking, no shared memory, no IPC. See `COMPLIANCE.md`.

## Tools

| Tool | Status | Purpose |
|------|--------|---------|
| `splitLayers` | working prototype (WP 5.4-spike) | Split wall-adjacent hex cells into geometrically graded boundary layers: ring-point insertion along columns, smoothed-normal or edge placement (`-placement`), offset-vector smoothing optimizer (`-optimizeSweeps`). Verified by the cavity smoke test in all ablation modes |
| `layerZoneExtract` | working (Stage A, 2026-06-18) | Emit protected cellZones (+ dynamicMeshDict fragment) for AMR from layer-cell reports — protects AMR-destructible extrusion layers (WP 5.3) |
| `layerAmrSurvival` | working (Stage A, 2026-06-18) | A/B AMR-survival test for boundary layers: refine+unrefine a sampled set, `-protect` freezes the layer cellZone from unrefinement (protectedCell_); decides AMR-native vs needs-protection (I4.2 guard) |
| `refinementFieldWriter` | planned (WP 7.2) | Error-indicator fields for dynamicRefineFvMesh |

## Building

Requires a sourced OpenFOAM environment (tested against OpenFOAM v2406+).

```sh
# WSL2 / any sourced OpenFOAM shell
wmake applications/splitLayers
```

CI builds run inside the `opencfd/openfoam-default` Docker image — see
`.github/workflows/build.yml`.

## Layout

```
applications/   one directory per utility (wmake-built)
etc/contracts/  versioned JSON schemas for spec files (semver, D7)
tests/smoke/    minimal cases exercising each tool end-to-end
```

## License

GPL-3.0-or-later. See `LICENSE`.
Vendored third-party code: `applications/splitLayers/json.hpp`
(nlohmann/json v3.11.3, MIT license — GPL-compatible).
OpenFOAM is a registered trademark of OpenCFD Ltd.
