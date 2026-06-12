# Compliance posture

This repository exists to keep the licensing boundary of the Zefra meshing
engine clean, auditable, and public (engine plan decisions D7, D10, attorney
review per D11).

## Rules

1. **Everything that links OpenFOAM/cfMesh lives here, under GPL-3.0, public
   from day one.** No OpenFOAM-linked binary is ever distributed from the
   proprietary engine repo.
2. **Process boundary only.** Tools communicate with the engine exclusively
   through versioned file contracts (`etc/contracts/`, JSON, semver,
   schema-validated on both sides). No linking against proprietary code, no
   shared memory, no sockets.
3. **No proprietary intelligence.** Tools execute fully-specified plans; they
   make no meshing decisions. If a change would move decision logic into a
   tool here, it belongs in the engine instead — redesign the contract.
4. **Contract changes are semver'd.** Breaking spec changes bump the major
   version; both sides validate the version field before execution.
5. **Third-party additions** to this repo must be GPL-compatible. The
   proprietary engine's own dependency allow-list is enforced separately
   (engine repo CI, decision D9).

## Review trail

- Architecture (process-boundary) opinion: attorney review scheduled
  Phase 2, Weeks 8–11 (D11). Opinion to be kept on file, referenced here.
- CGAL Alpha Wrapping commercial license (D1): separate executor
  (`wrap-cgal`), not part of this repository.
