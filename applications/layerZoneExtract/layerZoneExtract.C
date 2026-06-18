/*---------------------------------------------------------------------------*\
    layerZoneExtract  —  WP 5.3 AMR protection for extrusion boundary layers

    Extrusion (snappy addLayers) layers are AMR-DESTRUCTIBLE — unlike split-based
    layers, which are AMR-native (the cycle is perfectly reversible, I0.2). The
    internal-flow / manifold case ships on extrusion (the dual-executor gate), so
    its layer cells must be PROTECTED from AMR unrefinement or an AMR pass destroys
    them.

    This is a MECHANICAL executor (zero proprietary intelligence, D10): it reads
    the layer-cell set the proprietary executor wrote (meshcore.executors.
    snappy_layers.compute_layer_cells — the base-vs-layered nCells diff) and emits
    a protected cellZone (+ an optional dynamicMeshDict fragment) the AMR setup
    freezes. "Which cells are layers" is decided on the proprietary side; this tool
    only converts a set into a zone.

    DoD (the §8 headline): with the emitted zone protected, a refine/unrefine pass
    on the extrusion-layered manifold destroys ZERO layer cells (use
    refineUnrefineCycle as the deterministic check; a stock dynamicRefineFvMesh run
    is the real-pathway confirmation).

    Usage:  layerZoneExtract <cellSet> [-zone <name>] [-dynamicMeshDict]
    License: GPL-3.0-or-later (links OpenFOAM libraries).
\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "cellZone.H"
#include "OFstream.H"
#include "OSspecific.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Emit a protected cellZone (+ optional dynamicMeshDict fragment) from a "
        "boundary-layer cell set, to protect extrusion layers under AMR (WP 5.3)."
    );
    argList::addArgument("cellSet", "name of the layer-cell set (e.g. layerCells)");
    argList::addOption("zone", "name", "output cellZone name (default = set name)");
    argList::addBoolOption
    (
        "dynamicMeshDict",
        "also write system/dynamicMeshDict.layerProtect (a protection template)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word setName(args[1]);
    const word zoneName(args.getOrDefault<word>("zone", setName));

    cellSet layerCells(mesh, setName);
    const labelList cellLabels(layerCells.toc());
    Info<< "LZE|input|set|" << setName << "|cells|" << cellLabels.size() << nl;

    if (cellLabels.empty())
    {
        WarningInFunction
            << "cellSet '" << setName << "' is empty — no layer cells to protect "
            << "(snappy may have hit the first-layer/wall-cell floor)." << nl;
    }

    // Build (or replace) the cellZone from the set — the setsToZones idiom.
    // NB cellZoneMesh::set / the cellZone ctor are the version-sensitive lines;
    // verify against OpenFOAM v2512 on the first wmake.
    cellZoneMesh& zones = mesh.cellZones();
    label zid = zones.findZoneID(zoneName);
    if (zid < 0)
    {
        zid = zones.size();
        zones.setSize(zid + 1);
        zones.set
        (
            zid,
            new cellZone(zoneName, cellLabels, zid, zones)
        );
    }
    else
    {
        zones[zid] = cellLabels;   // cellZone::operator=(const labelUList&)
    }
    zones.clearAddressing();

    // Persist constant/polyMesh/cellZones. CRITICAL: a programmatically-added
    // cellZone is NOT written by mesh.write() unless the cellZoneMesh is marked
    // AUTO_WRITE — its default writeOpt is NO_WRITE when the mesh was read WITHOUT
    // zones, so mesh.write() silently skips it (the bug that made the zone never
    // reach disk while rc stayed 0). Set AUTO_WRITE, then write the zones directly.
    mesh.cellZones().writeOpt(IOobject::AUTO_WRITE);
    const bool ok = mesh.cellZones().write();
    Info<< "LZE|zone|" << zoneName << "|id|" << zid
        << "|cells|" << cellLabels.size() << "|written|" << ok << nl;

    if (args.found("dynamicMeshDict"))
    {
        mkDir(runTime.system());
        OFstream os(runTime.system()/"dynamicMeshDict.layerProtect");
        os  << "// layerZoneExtract — protect the '" << zoneName << "' layer cells\n"
            << "// from AMR unrefinement (hold them at their current refinement\n"
            << "// level / exclude from the unrefine candidate set). Merge into\n"
            << "// system/dynamicMeshDict; finalize exact keys on WSL (OpenFOAM\n"
            << "// v2512 dynamicRefineFvMesh).\n\n"
            << "dynamicFvMesh   dynamicRefineFvMesh;\n\n"
            << "dynamicRefineFvMeshCoeffs\n{\n"
            << "    // ... refine field / levels / intervals ...\n"
            << "    // freeze the boundary-layer cells:\n"
            << "    protectedZones ( " << zoneName << " );\n"
            << "}\n";
        Info<< "LZE|dynamicMeshDict|written|system/dynamicMeshDict.layerProtect" << nl;
    }

    Info<< "LZE|RESULT|zone|" << zoneName << "|cells|" << cellLabels.size()
        << "|ok|" << ok << nl;
    Info<< "End\n" << endl;
    return 0;
}
