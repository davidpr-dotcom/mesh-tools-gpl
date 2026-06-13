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

- Architecture (process-boundary) opinion: attorney review (D11), re-scoped per
  R8.3 to also cover GPL-oracle internal use + the clean-room protocol. Status,
  attorney questions, and the full third-party license register live in
  `Zefra/2026-06-13-trackB-legal-checkpoint.md`. Opinion to be kept on file,
  referenced here.
- CGAL Alpha Wrapping (D1 → **D1-rev**, 2026-06-12): the commercial license is
  **not** pursued. CGAL is used only as an INTERNAL, never-distributed black-box
  test ORACLE — clean-room (public API + the published paper only; the
  `Alpha_wrap` source is never read), kept outside this repo and the proprietary
  engine. The MVP wrapper is OpenVDB (`wrap-vdb`, MPL-2.0); the in-house
  `wrap-zefra` (Geogram/BSD) is the Track-W path to a proprietary primary.
