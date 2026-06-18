/*---------------------------------------------------------------------------*\
    layerAmrSurvival  —  A/B AMR-survival test for boundary-layer (BL) cells

    Question: under AMR, do the boundary-layer cells survive, and does PROTECTING
    them (layerZoneExtract's cellZone) guarantee it? This is the experiment that
    decides whether the dual-executor's layer meshes need protection or are
    AMR-native on their own (the I0.2 finding for split; measured equal for
    extrusion on the manifold).

    Method: refine a (sampled, hex-only) candidate set with hexRef8 — the SAME
    octree engine dynamicRefineFvMesh uses — then unrefine back. The candidate set
    is identical in both A/B modes; the ONLY difference is `-protect`, which mirrors
    dynamicRefineFvMesh's `protectedCell_`: protected cells are excluded from
    UNREFINEMENT (2:1 may still refine them — protectedCell_ blocks coarsening, not
    refinement). So:

      unprotected  blStuck = BL-zone cells NOT back at level 0 = reversibility
                             failures = the real "AMR destroyed the BL" signal
                             (manifold measured 0 → AMR-native).
      protected    blStuck = BL-zone cells that 2:1-refined and were KEPT (protected
                             from coarsening) = protection active, not destruction.

    Sampling caps the refined set so million-cell layer bands don't OOM.

    Usage:  layerAmrSurvival [-zone <name>] [-protect] [-sample <N>]
    License: GPL-3.0-or-later (links OpenFOAM libraries).
\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "hexRef8.H"
#include "hexMatcher.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "cellZone.H"
#include "bitSet.H"
#include "DynamicList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A/B AMR-survival test: refine a sampled set then unrefine; -protect freezes "
        "the layer cellZone from UNREFINEMENT (protectedCell_ semantics). Compare "
        "-protect 0 vs 1 to decide whether the BL needs protection."
    );
    argList::addOption("zone", "name", "protected layer cellZone (default: layerCells)");
    argList::addBoolOption("protect", "exclude the zone from unrefinement (protectedCell_)");
    argList::addOption("sample", "N", "cap the refined candidate count (default 50000)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word  zoneName = args.getOrDefault<word>("zone", "layerCells");
    const bool  protect  = args.found("protect");
    const label sampleCap = args.getOrDefault<label>("sample", 50000);

    const label zid = mesh.cellZones().findZoneID(zoneName);
    if (zid < 0)
    {
        FatalErrorInFunction
            << "cellZone '" << zoneName << "' not found — run layerZoneExtract first"
            << exit(FatalError);
    }
    const label nZone0 = mesh.cellZones()[zid].size();

    hexRef8 meshCutter(mesh);
    hexMatcher hexTest;

    // Candidate refine set: a deterministic strided sample of TRUE hexes (hexRef8
    // aborts on non-hexes). Identical in both A/B modes.
    const label nC0 = mesh.nCells();
    const label stride = (sampleCap > 0 && nC0 > sampleCap) ? (nC0 / sampleCap) : 1;
    DynamicList<label> cand(min(sampleCap, nC0));
    label nonHex = 0;
    for (label c = 0; c < nC0; c += stride)
    {
        if (hexTest.isA(mesh, c)) { cand.append(c); } else { ++nonHex; }
    }
    cand.shrink();
    Info<< "LAS|candidates|" << cand.size() << "|nonHexSkipped|" << nonHex
        << "|protect|" << label(protect) << "|zoneCells|" << nZone0
        << "|stride|" << stride << nl;

    const label nCells0 = mesh.nCells();
    const labelList toRefine = meshCutter.consistentRefinement(cand, true);

    // ---- REFINE ----
    {
        polyTopoChange meshMod(mesh);
        meshCutter.setRefinement(toRefine, meshMod);
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());
        meshCutter.updateMesh(map());
    }
    Info<< "LAS|afterRefine|cells|" << mesh.nCells() << nl;

    // ---- UNREFINE back ----  (protectedCell_: skip split points adjacent to the
    // protected zone, which auto-updated through the refine map)
    labelList splitPoints = meshCutter.getSplitPoints();
    if (protect)
    {
        bitSet prot(mesh.nCells(), false);
        for (const label c : mesh.cellZones()[zid]) { prot.set(c); }
        const labelListList& pCells = mesh.pointCells();
        DynamicList<label> keep(splitPoints.size());
        forAll(splitPoints, i)
        {
            const label pt = splitPoints[i];
            bool touchesProtected = false;
            for (const label c : pCells[pt])
            {
                if (prot.test(c)) { touchesProtected = true; break; }
            }
            if (!touchesProtected) { keep.append(pt); }
        }
        keep.shrink();
        Info<< "LAS|protect|splitPointsKept|" << keep.size()
            << "|ofTotal|" << splitPoints.size() << nl;
        splitPoints = keep;
    }
    {
        polyTopoChange meshMod(mesh);
        meshCutter.setUnrefinement(splitPoints, meshMod);
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());
        meshCutter.updateMesh(map());
    }
    Info<< "LAS|afterUnrefine|cells|" << mesh.nCells() << nl;

    // ---- verdict ----
    const labelList& cl = meshCutter.cellLevel();
    label blStuck = 0;
    for (const label c : mesh.cellZones()[zid]) { if (cl[c] != 0) { ++blStuck; } }
    label globalStuck = 0;
    forAll(cl, c) { if (cl[c] != 0) { ++globalStuck; } }

    Info<< "LAS|RESULT"
        << "|protect|"     << label(protect)
        << "|zoneCells|"   << mesh.cellZones()[zid].size()
        << "|candidates|"  << cand.size()
        << "|cellDelta|"   << (mesh.nCells() - nCells0)   // 0 = fully reversible
        << "|blStuck|"     << blStuck                     // unprotected: reversibility failures
        << "|globalStuck|" << globalStuck
        << nl;
    Info<< "End\n" << endl;
    return 0;
}
