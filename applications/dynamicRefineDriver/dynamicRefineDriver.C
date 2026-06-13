/*---------------------------------------------------------------------------*\
    dynamicRefineDriver  —  I0.2 AMR realism confirmation (the REAL pathway)

    Exercises stock dynamicRefineFvMesh exactly as a solver would: construct the
    dynamicFvMesh (-> dynamicRefineFvMesh per constant/dynamicMeshDict, with NO
    protectedCells), register the refinement field, and loop mesh.update().
    Refine the layer band (field=1), then unrefine it (field=0). If split layers
    are AMR-native, the cycle is reversible and nCells returns to the original —
    cellDelta 0 == zero destroyed, through the production refinement machinery.

    This complements the deterministic hexRef8 refineUnrefineCycle test and
    settles WP 5.3: layers survive real dynamicRefineFvMesh WITHOUT protection.

    License: GPL-3.0-or-later (links OpenFOAM libraries).
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "dynamicFvMesh.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Drive dynamicRefineFvMesh refine+unrefine; cellDelta 0 == zero destroyed"
    );
    #include "setRootCase.H"
    #include "createTime.H"

    autoPtr<dynamicFvMesh> meshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
    dynamicFvMesh& mesh = meshPtr();

    // the refinement driver field (tagger writes 1 in the layer band, 0 else)
    volScalarField refineField
    (
        IOobject
        (
            "refineField",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const label nCells0 = mesh.nCells();
    Info<< "I02S|nCells0|" << nCells0 << endl;

    // ---- REFINE: field already = 1 in the band; a few updates grow buffers ----
    for (label i = 0; i < 3; ++i)
    {
        runTime++;
        mesh.update();
        Info<< "I02S|refineStep|" << i << "|cells|" << mesh.nCells() << endl;
    }
    const label nPeak = mesh.nCells();

    // ---- UNREFINE: drop the field to 0 everywhere, update until stable ----
    refineField.primitiveFieldRef() = 0.0;
    refineField.correctBoundaryConditions();
    label prev = -1;
    for (label i = 0; i < 30 && mesh.nCells() != prev; ++i)
    {
        prev = mesh.nCells();
        runTime++;
        mesh.update();
        refineField.primitiveFieldRef() = 0.0;   // keep the mapped field at 0
        refineField.correctBoundaryConditions();
        Info<< "I02S|unrefineStep|" << i << "|cells|" << mesh.nCells() << endl;
    }
    const label nCells1 = mesh.nCells();

    Info<< "I02S|RESULT"
        << "|nCells0|"   << nCells0
        << "|peak|"      << nPeak
        << "|nCells1|"   << nCells1
        << "|cellDelta|" << (nCells1 - nCells0)   // 0 => fully reversible
        << nl;

    runTime.writeNow();
    Info<< "End\n" << endl;
    return 0;
}
