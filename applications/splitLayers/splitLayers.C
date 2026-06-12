/*---------------------------------------------------------------------------*\
    splitLayers — split-based boundary-layer generation (WP 5.4-spike)

    Splits wall-adjacent cells into n geometrically graded layers toward the
    wall, with an optional node-position optimization pass.

    STATUS: prototype stub. Current behavior: parse arguments, read the
    layer-spec file, load the mesh, identify wall-adjacent cells for each
    patch named in the spec, report, exit. The split mechanism lands next
    (spike Week 8), optimization pass after (Week 9).

    Contract: consumes a layer-spec file (etc/contracts/layer-spec
    schema, semver-checked). Emits structured diagnostics on stdout with
    the `SPLITLAYERS|` prefix for machine parsing.

    This tool contains no planning logic: layer counts, heights and
    expansion ratios arrive fully specified. See COMPLIANCE.md.

    License: GPL-3.0-or-later (links OpenFOAM libraries).
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Split wall-adjacent cells into graded boundary layers "
        "(prototype, WP 5.4-spike)"
    );

    argList::addOption
    (
        "spec",
        "file",
        "Layer-spec file (JSON, contract version layer-spec/0.x)"
    );

    argList::addBoolOption
    (
        "dryRun",
        "Identify and report wall-adjacent cells only; modify nothing"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Create mesh for time = " << runTime.timeName() << nl << endl;

    fvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // ---- Spec file (presence check only, parsing lands with the split) ----
    fileName specFile;
    if (args.readIfPresent("spec", specFile))
    {
        IFstream is(specFile);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot open spec file " << specFile << exit(FatalError);
        }
        Info<< "SPLITLAYERS|spec|" << specFile << nl;
        // TODO(WP5.4-spike): parse JSON (vendored MIT single-header parser),
        // validate semver against supported contract range, extract per-patch
        // layer specs {nLayers, firstLayerHeight, expansionRatio}.
    }
    else
    {
        Info<< "SPLITLAYERS|spec|none (reporting all wall patches)" << nl;
    }

    const bool dryRun = args.found("dryRun");

    // ---- Wall-adjacent cell identification ----
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() || isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        // Owner cells of boundary faces = the first layer front
        const labelUList& faceCells = pp.faceCells();

        labelHashSet adjacentCells(faceCells.size());
        forAll(faceCells, i)
        {
            adjacentCells.insert(faceCells[i]);
        }

        Info<< "SPLITLAYERS|patch|" << pp.name()
            << "|nFaces|" << pp.size()
            << "|nAdjacentCells|" << adjacentCells.size()
            << nl;

        // TODO(WP5.4-spike): per-cell wall-normal direction from owner-face
        // normal, smoothed over the layer front; directional split x n with
        // geometric spacing toward the wall; then node-position optimization
        // (condition-number objective, wall nodes fixed).
    }

    if (dryRun)
    {
        Info<< "SPLITLAYERS|result|dry-run-ok" << nl;
    }
    else
    {
        Info<< "SPLITLAYERS|result|no-op (split mechanism not yet implemented)"
            << nl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
