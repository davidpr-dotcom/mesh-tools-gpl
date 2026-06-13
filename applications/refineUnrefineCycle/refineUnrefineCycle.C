/*---------------------------------------------------------------------------*\
    refineUnrefineCycle  —  I0.2 AMR structural test (gate criterion 3)

    Refine a cellSet (the layer band) then unrefine it back, using hexRef8 —
    the SAME octree engine dynamicRefineFvMesh uses at runtime. If split-based
    layers are ordinary hexes (the spike's structural promise), the cycle is
    perfectly reversible: nCells returns, every cell drops back to level 0, and
    ZERO layer cells are destroyed. This is the deterministic CI metric; a stock
    dynamicRefineFvMesh run confirms the same through the real AMR pathway.

    Usage:  refineUnrefineCycle <cellSetName>
    License: GPL-3.0-or-later (links OpenFOAM libraries).
\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "hexRef8.H"
#include "hexMatcher.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "cellSet.H"
#include "DynamicList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Refine + unrefine a cellSet (layer band) and count destroyed cells"
    );
    argList::addArgument("cellSet", "name of the cellSet to cycle");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word setName = args[1];

    const label  nCells0 = mesh.nCells();
    const label  nFaces0 = mesh.nFaces();
    const scalar vol0    = gSum(mesh.cellVolumes());

    // hexRef8 reads constant/polyMesh/{cellLevel,pointLevel} if present
    // (the tagger writes all-zeros), else initialises to level 0.
    hexRef8 meshCutter(mesh);

    cellSet bandSet(mesh, setName);
    const labelList raw(bandSet.toc());

    // hexRef8 only refines TRUE hexahedra. The tagger's cfc==6 heuristic lets
    // through 6-face non-hexes (e.g. concave-corner cells on the bracket),
    // which abort hexRef8 -> filter with the real hex matcher here.
    hexMatcher hexTest;
    DynamicList<label> keep(raw.size());
    forAll(raw, i)
    {
        if (hexTest.isA(mesh, raw[i])) { keep.append(raw[i]); }
    }
    keep.shrink();
    const labelList candidate(keep);
    Info<< "I02|refineCandidates|" << candidate.size()
        << "|nonHexDropped|" << (raw.size() - candidate.size()) << nl;

    // 2:1 consistency (only refinable hexes are pulled in; non-hex base cells
    // are level-locked and never forced)
    labelList cellsToRefine = meshCutter.consistentRefinement(candidate, true);
    Info<< "I02|consistentRefine|" << cellsToRefine.size() << nl;

    // ---- REFINE ----
    {
        polyTopoChange meshMod(mesh);
        meshCutter.setRefinement(cellsToRefine, meshMod);
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());
        meshCutter.updateMesh(map());
    }
    Info<< "I02|afterRefine|cells|" << mesh.nCells() << nl;

    // ---- UNREFINE back (all eligible split points) ----
    labelList splitPoints = meshCutter.getSplitPoints();
    Info<< "I02|splitPoints|" << splitPoints.size() << nl;
    {
        polyTopoChange meshMod(mesh);
        meshCutter.setUnrefinement(splitPoints, meshMod);
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
        mesh.updateMesh(map());
        meshCutter.updateMesh(map());
    }
    Info<< "I02|afterUnrefine|cells|" << mesh.nCells() << nl;

    // ---- verdict ----
    const label  nCells1 = mesh.nCells();
    const label  nFaces1 = mesh.nFaces();
    const scalar vol1    = gSum(mesh.cellVolumes());
    const labelList& cl  = meshCutter.cellLevel();
    label nStuck = 0;
    forAll(cl, c) { if (cl[c] != 0) { ++nStuck; } }

    Info<< "I02|RESULT"
        << "|nCells0|"      << nCells0
        << "|nCells1|"      << nCells1
        << "|cellDelta|"    << (nCells1 - nCells0)        // must be 0
        << "|cellsStuckRefined|" << nStuck                // must be 0 (= destroyed)
        << "|nFacesDelta|"  << (nFaces1 - nFaces0)        // must be 0
        << "|volErr|"       << mag(vol1 - vol0)/max(vol0, VSMALL)
        << nl;

    // write the cycled mesh so the Python comparator can confirm checkMesh
    // metric invariance (the DoD 'harness metrics unchanged')
    mesh.write();
    Info<< "End\n" << endl;
    return 0;
}
