/*---------------------------------------------------------------------------*\
    splitLayers — split-based boundary-layer generation (WP 5.4-spike)

    Multi-row columns (see docs/multi-row-columns.md): a column is a CHAIN
    of hex cells walked from the wall face through successive opposing
    faces, so the layer stack is no longer bounded by the (possibly
    trimmed) wall-cell height. Ring points are placed along the chain's
    four corner RAILS at arc-length positions; the chain's interface faces
    are absorbed; side faces are partitioned by the rings assigned to
    their segment; rail joints survive as polygon vertices for lateral
    neighbours.

    Consistency (v1, strictness over cleverness): rings are a property of
    the rail (first column wins); columns whose shared rails disagree on
    scale/segmentation, whose own rails disagree on ring->segment
    assignment, or whose side faces are claimed by a crossing chain are
    SKIPPED and reported. All conflicts are detected before any point is
    created.

    Stages: spec (layer-spec JSON 0.x) -> identification -> conflict
    pre-passes -> ring creation -> topology (polyTopoChange) -> optional
    offset-vector smoothing (-optimizeSweeps, wall points fixed).
    `-placement normal` applies to single-cell columns only.

    Diagnostics on stdout with the `SPLITLAYERS|` prefix.
    License: GPL-3.0-or-later (links OpenFOAM libraries).
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
#include "EdgeMap.H"
#include "edgeList.H"
#include "unitConversion.H"
#include "OFstream.H"

#include <algorithm>
#include <fstream>
#include <string>
#include "json.hpp"

using namespace Foam;
using nlohmann::json;

// * * * * * * * * * * * * * * * * * spec  * * * * * * * * * * * * * * * * //

struct PatchSpec
{
    word patchName;
    label nLayers;
    scalar firstHeight;
    scalar expansion;
};

static List<PatchSpec> readSpec(const fileName& specFile)
{
    std::ifstream is(specFile.c_str());
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open spec file " << specFile << exit(FatalError);
    }
    json j = json::parse(is);

    if (j.value("contract", "") != std::string("layer-spec"))
    {
        FatalErrorInFunction
            << "spec contract is not 'layer-spec'" << exit(FatalError);
    }
    const std::string ver = j.value("version", "");
    if (ver.rfind("0.", 0) != 0)
    {
        FatalErrorInFunction
            << "unsupported layer-spec major version: " << ver.c_str()
            << " (tool supports 0.x)" << exit(FatalError);
    }
    Info<< "SPLITLAYERS|spec|version|" << ver.c_str() << nl;

    List<PatchSpec> specs(j["patches"].size());
    label i = 0;
    for (const auto& p : j["patches"])
    {
        specs[i].patchName = word(p["name"].get<std::string>());
        specs[i].nLayers = p["nLayers"].get<int>();
        specs[i].firstHeight = p["firstLayerHeight"].get<double>();
        specs[i].expansion = p["expansionRatio"].get<double>();
        ++i;
    }
    return specs;
}

// * * * * * * * * * *  smoothed inward point normals  * * * * * * * * * * //

static vectorField smoothedInwardNormals
(
    const polyPatch& pp,
    const label nSmooth,
    const scalar featureAngleDeg,
    label& nFeaturePoints
)
{
    const vectorField faceN(pp.faceNormals());     // outward of domain
    const labelListList& pointFaces = pp.pointFaces();

    vectorField n(pp.nPoints(), Zero);
    forAll(pointFaces, pi)
    {
        for (const label fi : pointFaces[pi])
        {
            n[pi] -= faceN[fi] * pp.magFaceAreas()[fi];   // inward
        }
    }
    n /= max(mag(n), SMALL);

    const scalar cosFeat = Foam::cos(degToRad(featureAngleDeg));
    nFeaturePoints = 0;
    forAll(pointFaces, pi)
    {
        const labelList& pf = pointFaces[pi];
        bool feat = false;
        forAll(pf, a)
        {
            for (label b = a + 1; b < pf.size(); ++b)
            {
                if ((faceN[pf[a]] & faceN[pf[b]]) < cosFeat) { feat = true; }
            }
        }
        if (feat) { ++nFeaturePoints; }
    }

    const edgeList& edges = pp.edges();
    for (label it = 0; it < nSmooth; ++it)
    {
        vectorField acc(n.size(), Zero);
        scalarField cnt(n.size(), 0.0);
        for (const edge& e : edges)
        {
            acc[e[0]] += n[e[1]];  cnt[e[0]] += 1.0;
            acc[e[1]] += n[e[0]];  cnt[e[1]] += 1.0;
        }
        forAll(n, pi)
        {
            if (cnt[pi] > 0)
            {
                n[pi] = 0.5 * n[pi] + 0.5 * acc[pi] / cnt[pi];
            }
        }
        n /= max(mag(n), SMALL);
    }
    return n;
}

// * * * * * * * * * * * * * column structures * * * * * * * * * * * * * * //

struct Column
{
    label patchFacei = -1;
    label meshFacei = -1;       // wall face
    label si = -1;              // spec index
    labelList chain;            // cells, wall first
    labelList ifaces;           // interfaces between chain cells (m-1)
    label topFace = -1;         // opposing face of the last cell
    FixedList<labelList, 4> rails;   // point paths, size m+1 each
    scalar scale = 1.0;
    bool skip = false;
    const char* reason = nullptr;
};

struct RailData
{
    labelList joints;           // FULL point path — the rail's identity.
                                // Keying by first edge alone is ambiguous:
                                // at corner bisectors two columns share
                                // their first rail edge and then DIVERGE
                                // (2026-06-12 forensics: strips built with
                                // rings of the other path's rings).
    labelList rings;            // n new point labels, wall-side first
    labelList assign;           // n: segment index of each ring
    label nSeg = 0;
    scalar scale = 1.0;
    label si = -1;
};

int main(int argc, char *argv[])
{
    argList::addNote("Split-based boundary layers (multi-row columns)");
    argList::addOption("spec", "file", "layer-spec JSON (contract 0.x)");
    argList::addOption("placement", "edge|normal",
                       "single-cell-column placement (default: normal)");
    argList::addOption("optimizeSweeps", "N",
                       "tangential smoothing sweeps (default 10, 0=off)");
    argList::addOption("featureAngle", "deg",
                       "feature-point detection angle (default 45)");
    argList::addBoolOption("dryRun", "identify and report columns only");
    argList::addBoolOption("debugDump",
        "write constant/splitLayersDebugMap: finalFace origFace tag column");
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Create mesh for time = " << runTime.timeName() << nl << endl;
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    fileName specFile;
    if (!args.readIfPresent("spec", specFile))
    {
        FatalErrorInFunction << "-spec is required" << exit(FatalError);
    }
    const List<PatchSpec> specs = readSpec(specFile);
    const word placement = args.getOrDefault<word>("placement", "normal");
    const label nOptSweeps = args.getOrDefault<label>("optimizeSweeps", 10);
    const scalar featAngle = args.getOrDefault<scalar>("featureAngle", 45.0);
    const bool dryRun = args.found("dryRun");
    const bool debugDump = args.found("debugDump");

    // forensics: original-face -> (construct tag, column id)
    Map<word> faceTagOf;
    Map<label> faceColOf;
    auto recordFace = [&](const label origFace, const char* tag,
                          const label ci2)
    {
        if (!debugDump) { return; }
        if (faceTagOf.found(origFace))
        {
            // multiple constructs touching one original face is itself
            // suspicious — keep both tags
            faceTagOf[origFace] = faceTagOf[origFace] + "+" + tag;
        }
        else
        {
            faceTagOf.insert(origFace, word(tag));
            faceColOf.insert(origFace, ci2);
        }
    };

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelListList& cellEdges = mesh.cellEdges();
    const edgeList& meshEdges = mesh.edges();
    const pointField& pts = mesh.points();

    // heights
    List<scalarList> cumHeights(specs.size());
    scalarList needed(specs.size());
    forAll(specs, si)
    {
        cumHeights[si].setSize(specs[si].nLayers);
        scalar t = 0, h = specs[si].firstHeight;
        forAll(cumHeights[si], k)
        {
            t += h;
            cumHeights[si][k] = t;
            h *= specs[si].expansion;
        }
        // healthy remainder: continue the progression one more step
        needed[si] = t + specs[si].firstHeight
            * Foam::pow(specs[si].expansion, specs[si].nLayers);
    }

    labelList patchFaceCount(mesh.nCells(), 0);
    forAll(specs, si)
    {
        const label patchi = patches.findPatchID(specs[si].patchName);
        if (patchi < 0)
        {
            FatalErrorInFunction
                << "patch '" << specs[si].patchName << "' not found"
                << exit(FatalError);
        }
        for (const label c : patches[patchi].faceCells())
        {
            patchFaceCount[c]++;
        }
    }

    // helper: is this cell a clean hex (6 quad faces)?
    auto isHex = [&](const label celli) -> bool
    {
        const cell& c = cells[celli];
        if (c.size() != 6) { return false; }
        for (const label fi : c)
        {
            if (faces[fi].size() != 4) { return false; }
        }
        return true;
    };

    // helper: corner continuation pt -> next point on target face via a
    // cell edge of celli
    auto railStep = [&](const label celli, const label fromPt,
                        const face& target) -> label
    {
        for (const label ei : cellEdges[celli])
        {
            const edge& e = meshEdges[ei];
            const label other = e.otherVertex(fromPt);
            if (other >= 0 && target.found(other)) { return other; }
        }
        return -1;
    };

    // ---- PASS 1: identification ----------------------------------------
    DynamicList<Column> columns;
    labelList claimed(mesh.nCells(), -1);
    List<label> nCorner(specs.size(), 0), nNonHex(specs.size(), 0);

    forAll(specs, si)
    {
        const polyPatch& pp = patches[patches.findPatchID(specs[si].patchName)];

        forAll(pp, pfi)
        {
            const label fG = pp.start() + pfi;
            const label own = mesh.faceOwner()[fG];

            if (patchFaceCount[own] > 1) { ++nCorner[si]; continue; }
            if (!isHex(own) || claimed[own] >= 0)
            {
                ++nNonHex[si];
                continue;
            }

            Column col;
            col.patchFacei = pfi;
            col.meshFacei = fG;
            col.si = si;

            const face& wf = faces[fG];
            label iface = cells[own].opposingFaceLabel(fG, faces);
            if (iface < 0 || faces[iface].size() != 4)
            {
                ++nNonHex[si];
                continue;
            }

            bool ok = true;
            forAll(wf, c)
            {
                col.rails[c].append(wf[c]);
                const label q = railStep(own, wf[c], faces[iface]);
                if (q < 0) { ok = false; break; }
                col.rails[c].append(q);
            }
            if (!ok) { ++nNonHex[si]; continue; }

            col.chain.append(own);
            claimed[own] = columns.size();

            // extend while the shortest rail is shorter than needed
            auto minRailLen = [&]() -> scalar
            {
                scalar L = GREAT;
                forAll(col.rails, c)
                {
                    scalar s = 0;
                    const labelList& r = col.rails[c];
                    for (label j = 1; j < r.size(); ++j)
                    {
                        s += mag(pts[r[j]] - pts[r[j-1]]);
                    }
                    L = min(L, s);
                }
                return L;
            };

            while (minRailLen() < needed[si] && mesh.isInternalFace(iface))
            {
                const label cOwn = mesh.faceOwner()[iface];
                const label cNei = mesh.faceNeighbour()[iface];
                const label next = (cOwn == col.chain.last()) ? cNei : cOwn;
                if (next < 0) { break; }
                if (claimed[next] >= 0 || patchFaceCount[next] > 0
                 || !isHex(next))
                {
                    break;
                }
                const label nextIface =
                    cells[next].opposingFaceLabel(iface, faces);
                if (nextIface < 0 || faces[nextIface].size() != 4) { break; }

                bool stepOk = true;
                FixedList<label, 4> qs;
                forAll(col.rails, c)
                {
                    qs[c] = railStep(next, col.rails[c].last(),
                                     faces[nextIface]);
                    if (qs[c] < 0) { stepOk = false; break; }
                }
                if (!stepOk) { break; }

                col.ifaces.append(iface);
                col.chain.append(next);
                claimed[next] = columns.size();
                forAll(col.rails, c) { col.rails[c].append(qs[c]); }
                iface = nextIface;
            }
            col.topFace = iface;
            col.scale = min(scalar(1), minRailLen() / needed[si]);
            columns.append(col);
        }
    }

    // ---- PASS 2: conflict pre-passes ------------------------------------
    // Rail identity = FIRST SEGMENT EDGE (base point + direction). Keying
    // by base point alone is wrong at reentrant corner lines, where one
    // wall point carries rails of two perpendicular columns (found
    // 2026-06-12: teleported-vertex/bowtie faces along the bracket corner).
    EdgeMap<RailData> railsOf;
    auto railKey = [&](const Column& cl, const label c) -> edge
    {
        return edge(cl.rails[c][0], cl.rails[c][1]);
    };
    List<label> nConflict(specs.size(), 0), nWarp(specs.size(), 0),
                nCross(specs.size(), 0);

    // rail assignment for a column's rail c under the column's scale
    auto computeAssign = [&](const Column& col, const label c) -> labelList
    {
        const labelList& r = col.rails[c];
        scalarList arc(r.size(), Zero);
        for (label j = 1; j < r.size(); ++j)
        {
            arc[j] = arc[j-1] + mag(pts[r[j]] - pts[r[j-1]]);
        }
        const scalarList& cum = cumHeights[col.si];
        labelList assign(cum.size());
        forAll(cum, k)
        {
            const scalar h = col.scale * cum[k];
            label seg = 0;
            while (seg + 2 < r.size() && h > arc[seg + 1]) { ++seg; }
            assign[k] = seg;
        }
        return assign;
    };

    Map<label> sideFaceCol;   // side face -> column id (cross detection)
    Map<label> topFaceCol;    // top face -> column id (head-on detection)
    Map<label> ifaceCol;      // interface face -> column id (cross detection)

    forAll(columns, ci)
    {
        Column& col = columns[ci];
        if (col.skip) { continue; }
        const label m = col.chain.size();

        // (a) shared-rail consistency: chain length, scale AND the ring->
        // segment assignment table must match what a neighbouring column
        // already recorded on a shared rail. Assignment mismatch was the
        // bracket quality defect (2026-06-12 forensics: all 89 bad faces
        // were strips of columns whose shared rails were subdivided under
        // a neighbour's table while their own rails used their own).
        const labelList a0 = computeAssign(col, 0);
        bool conflict = false;
        forAll(col.rails, c)
        {
            const auto it = railsOf.cfind(railKey(col, c));
            if (it.good()
             && (it.val().nSeg != m
              || mag(it.val().scale - col.scale) > 1e-6 * col.scale + VSMALL
              || it.val().assign != a0
              || it.val().joints != col.rails[c]))   // FULL-PATH identity
            {
                conflict = true;
            }
        }
        if (conflict)
        {
            col.skip = true; col.reason = "conflict";
            ++nConflict[col.si];
            continue;
        }

        // (a2) degenerate chains: ANY duplicate point among the four rails
        // (converging rails in trimmed cells) breaks rail-pair
        // identification and point substitution — strips then borrow ring
        // points from a rail that is not on the face (2026-06-12
        // forensics: bowtie/teleported-vertex strips at mid-chain
        // segments). Skip strictly.
        {
            labelHashSet railPts(4 * (m + 1));
            bool degen = false;
            forAll(col.rails, c)
            {
                for (const label p : col.rails[c])
                {
                    if (!railPts.insert(p)) { degen = true; break; }
                }
                if (degen) { break; }
            }
            if (degen)
            {
                col.skip = true; col.reason = "warp";
                ++nWarp[col.si];
                continue;
            }
        }

        // (b) warp: all four rails must agree on ring->segment assignment
        bool warp = false;
        for (label c = 1; c < 4 && !warp; ++c)
        {
            if (computeAssign(col, c) != a0) { warp = true; }
        }
        if (warp)
        {
            col.skip = true; col.reason = "warp";
            ++nWarp[col.si];
            continue;
        }

        // (c) crossing chains + face classification. Every face of every
        // chain cell must be wall, interface, top, or a clean side face
        // (two rails' (i,i+1) point pairs). Side faces claimed by another
        // column must share rails (parallel chains); top faces meeting
        // another column are only allowed HEAD-ON (top-to-top). Anything
        // else: skip this column — a silently unhandled face leaves a
        // deleted-owner error in polyTopoChange.
        bool cross = false;
        bool unclass = false;
        forAll(col.chain, i)
        {
            for (const label fi : cells[col.chain[i]])
            {
                if (fi == col.meshFacei) { continue; }
                if (col.ifaces.found(fi))
                {
                    // our interface must not be anyone's side/top face
                    // (perpendicular chains: B strips A's interface ->
                    // A's remainder never closes — found 2026-06-12)
                    const auto its = sideFaceCol.cfind(fi);
                    const auto itt = topFaceCol.cfind(fi);
                    if ((its.good() && its.val() != ci)
                     || (itt.good() && itt.val() != ci))
                    {
                        cross = true;
                        break;
                    }
                    continue;
                }
                if (fi == col.topFace)
                {
                    // perpendicular meeting: our top is someone's side
                    // or interface
                    const auto its = sideFaceCol.cfind(fi);
                    const auto iti = ifaceCol.cfind(fi);
                    if ((its.good() && its.val() != ci)
                     || (iti.good() && iti.val() != ci))
                    {
                        cross = true;
                    }
                    continue;
                }

                // must be a clean side face
                label nr = 0;
                forAll(col.rails, c)
                {
                    if (faces[fi].found(col.rails[c][i])
                     && faces[fi].found(col.rails[c][i + 1]))
                    {
                        ++nr;
                    }
                }
                if (nr != 2) { unclass = true; break; }

                const auto it = sideFaceCol.cfind(fi);
                if (it.good() && it.val() != ci)
                {
                    const Column& other = columns[it.val()];
                    bool parallel = false;
                    forAll(col.rails, c)
                    {
                        if (faces[fi].found(col.rails[c][i]))
                        {
                            forAll(other.rails, oc)
                            {
                                if (other.rails[oc].found(col.rails[c][i]))
                                {
                                    parallel = true;
                                }
                            }
                        }
                    }
                    if (!parallel) { cross = true; }
                }
                // our side face is someone's TOP or INTERFACE face
                const auto itt = topFaceCol.cfind(fi);
                if (itt.good() && itt.val() != ci) { cross = true; }
                const auto iti = ifaceCol.cfind(fi);
                if (iti.good() && iti.val() != ci) { cross = true; }
            }
            if (cross || unclass) { break; }
        }
        if (unclass)
        {
            col.skip = true; col.reason = "warp";
            ++nWarp[col.si];
            continue;
        }
        if (cross)
        {
            col.skip = true; col.reason = "cross";
            ++nCross[col.si];
            continue;
        }
        forAll(col.chain, i)
        {
            for (const label fi : cells[col.chain[i]])
            {
                if (fi != col.meshFacei && fi != col.topFace
                 && !col.ifaces.found(fi))
                {
                    sideFaceCol.insert(fi, ci);
                }
            }
        }
        topFaceCol.insert(col.topFace, ci);
        for (const label fi : col.ifaces) { ifaceCol.insert(fi, ci); }

        // record rails (tentative ring labels filled in PASS 3)
        forAll(col.rails, c)
        {
            if (!railsOf.found(railKey(col, c)))
            {
                RailData rd;
                rd.joints = col.rails[c];
                rd.nSeg = m;
                rd.scale = col.scale;
                rd.si = col.si;
                rd.assign = a0;
                railsOf.insert(railKey(col, c), rd);
            }
        }
    }

    // report
    {
        List<label> nCols(specs.size(), 0), nMulti(specs.size(), 0);
        forAll(columns, ci)
        {
            if (!columns[ci].skip)
            {
                ++nCols[columns[ci].si];
                if (columns[ci].chain.size() > 1)
                {
                    ++nMulti[columns[ci].si];
                }
            }
        }
        forAll(specs, si)
        {
            Info<< "SPLITLAYERS|patch|" << specs[si].patchName
                << "|nFaces|"
                << patches[patches.findPatchID(specs[si].patchName)].size()
                << "|columns|" << nCols[si]
                << "|multiRow|" << nMulti[si]
                << "|skippedNonHex|" << nNonHex[si]
                << "|skippedCorner|" << nCorner[si]
                << "|skippedConflict|" << nConflict[si]
                << "|skippedWarp|" << nWarp[si]
                << "|skippedCross|" << nCross[si]
                << nl;
        }
    }

    if (dryRun)
    {
        Info<< "SPLITLAYERS|result|dry-run-ok" << nl << "End\n" << endl;
        return 0;
    }

    // smoothed normals per spec patch (single-cell-column placement)
    List<vectorField> patchNormals(specs.size());
    forAll(specs, si)
    {
        const polyPatch& pp = patches[patches.findPatchID(specs[si].patchName)];
        label nFeat = 0;
        patchNormals[si] = smoothedInwardNormals(pp, 5, featAngle, nFeat);
        Info<< "SPLITLAYERS|normals|" << specs[si].patchName
            << "|featurePoints|" << nFeat << nl;
    }

    // ---- PASS 3: ring creation ------------------------------------------
    polyTopoChange meshMod(mesh);
    EdgeMap<labelList> segRings;        // segment edge -> ordered ring pts
    EdgeMap<label> edgeClaim;           // segment edge -> owning rail base
    Map<point> newPtPos;                // topo point label -> coordinates

    label nEdgeConflict = 0;
    forAll(columns, ci)
    {
        Column& col0 = columns[ci];
        if (col0.skip) { continue; }

        // Pre-check: a rail not yet created must not reuse a segment edge
        // already claimed by a DIFFERENT rail (rails from different wall
        // points can merge onto the same mesh edges at reentrant corners
        // — found 2026-06-12: open-cell line along the bracket corner).
        bool edgeConflict = false;
        forAll(col0.rails, c)
        {
            if (railsOf[railKey(col0, c)].rings.size()) { continue; }
            const label base = col0.rails[c][0];
            const labelList& r = col0.rails[c];
            for (label j = 0; j + 1 < r.size() && !edgeConflict; ++j)
            {
                const auto it = edgeClaim.cfind(edge(r[j], r[j+1]));
                if (it.good() && it.val() != base) { edgeConflict = true; }
            }
            if (edgeConflict) { break; }
        }
        if (edgeConflict)
        {
            col0.skip = true;
            col0.reason = "conflict";
            ++nEdgeConflict;
            continue;
        }
        forAll(col0.rails, c)
        {
            const labelList& r = col0.rails[c];
            for (label j = 0; j + 1 < r.size(); ++j)
            {
                edgeClaim.insert(edge(r[j], r[j+1]), col0.rails[c][0]);
            }
        }

        const Column& col = columns[ci];
        const polyPatch& pp =
            patches[patches.findPatchID(specs[col.si].patchName)];
        const scalar cosGuard = Foam::cos(degToRad(60.0));

        forAll(col.rails, c)
        {
            const label base = col.rails[c][0];
            RailData& rd = railsOf[railKey(col, c)];
            if (rd.rings.size()) { continue; }   // already created

            const labelList& r = col.rails[c];
            scalarList arc(r.size(), Zero);
            for (label j = 1; j < r.size(); ++j)
            {
                arc[j] = arc[j-1] + mag(pts[r[j]] - pts[r[j-1]]);
            }

            // single-cell columns may use smoothed-normal placement
            vector dir = Zero;
            bool useDir = false;
            if (r.size() == 2 && placement == "normal")
            {
                const vector edgeDir =
                    (pts[r[1]] - pts[r[0]]) / max(arc.last(), VSMALL);
                const label loc = pp.whichPoint(base);
                if (loc >= 0)
                {
                    const vector& ns = patchNormals[col.si][loc];
                    if ((ns & edgeDir) > cosGuard)
                    {
                        dir = ns;
                        useDir = true;
                    }
                }
            }

            rd.rings.setSize(specs[col.si].nLayers);
            forAll(rd.rings, k)
            {
                const scalar h = rd.scale * cumHeights[col.si][k];
                point rp;
                if (useDir)
                {
                    rp = pts[r[0]] + dir * h;
                }
                else
                {
                    const label seg = rd.assign[k];
                    const scalar f =
                        (h - arc[seg])
                      / max(arc[seg + 1] - arc[seg], VSMALL);
                    rp = pts[r[seg]] + f * (pts[r[seg + 1]] - pts[r[seg]]);
                }
                rd.rings[k] = meshMod.addPoint(rp, base, -1, true);
                newPtPos.insert(rd.rings[k], rp);
            }

            // register rings per segment edge for the lateral pass
            forAll(r, j)
            {
                if (j + 1 >= r.size()) { break; }
                DynamicList<label> inSeg;
                forAll(rd.assign, k)
                {
                    if (rd.assign[k] == j) { inSeg.append(rd.rings[k]); }
                }
                if (inSeg.size())
                {
                    segRings.insert(edge(r[j], r[j+1]), labelList(inSeg));
                }
            }
        }
    }

    if (nEdgeConflict)
    {
        Info<< "SPLITLAYERS|lateSkips|edgeConflict|" << nEdgeConflict << nl;
    }

    // ---- PASS 4: topology -------------------------------------------------
    Map<labelList> stackOf;     // column id -> n+1 new cell labels
    Map<label> cellCol;         // chain cell -> column id
    forAll(columns, ci)
    {
        const Column& col = columns[ci];
        if (col.skip) { continue; }
        const label n = specs[col.si].nLayers;
        labelList stack(n + 1);
        forAll(stack, j)
        {
            stack[j] = meshMod.addCell(-1, -1, -1, col.chain[0], -1);
        }
        for (const label cc : col.chain)
        {
            meshMod.removeCell(cc, -1);
            cellCol.insert(cc, ci);
        }
        stackOf.insert(ci, stack);
    }

    bitSet faceDone(mesh.nFaces(), false);

    // Insert ring points of FOREIGN rails into a face polygon: at reentrant
    // geometry edges a split column's side/top face can contain another
    // column's rail segment as one of its edges (found 2026-06-12: open-
    // cell line along the bracket's concave corner). Same walk as the
    // lateral pass; own-rail edges never match segRings (their rings are
    // already face boundaries, and ringless segments have no entry).
    auto withForeignRings = [&](const face& f0) -> face
    {
        DynamicList<label> nf(f0.size() + 16);
        forAll(f0, k)
        {
            nf.append(f0[k]);
            const edge fe(f0[k], f0.nextLabel(k));
            const auto it = segRings.cfind(fe);
            if (it.good())
            {
                const labelList& ring = it.val();
                if (f0[k] == it.key()[0])
                {
                    forAll(ring, j) { nf.append(ring[j]); }
                }
                else
                {
                    forAllReverse(ring, j) { nf.append(ring[j]); }
                }
            }
        }
        return face(labelList(nf));
    };

    // construction-time sanity: every face we build must stay within the
    // size neighbourhood of the original face it derives from. Violations
    // identify the construct directly (debug aid; cheap enough to keep on).
    label nInsane = 0;
    auto saneCheck = [&](const face& nf, const label origFace,
                         const label ci2, const char* tag)
    {
        if (nInsane >= 20) { return; }
        const face& f0 = faces[origFace];
        scalar d0 = 0;
        forAll(f0, a)
        {
            for (label b2 = a + 1; b2 < f0.size(); ++b2)
            {
                d0 = max(d0, magSqr(pts[f0[a]] - pts[f0[b2]]));
            }
        }
        auto coordOf = [&](const label p) -> point
        {
            return (p < pts.size() && !newPtPos.found(p))
                 ? pts[p] : newPtPos[p];
        };
        scalar d1 = 0;
        forAll(nf, a)
        {
            for (label b2 = a + 1; b2 < nf.size(); ++b2)
            {
                d1 = max(d1, magSqr(coordOf(nf[a]) - coordOf(nf[b2])));
            }
        }
        if (d1 > 9.0 * d0 + SMALL)
        {
            ++nInsane;
            Info<< "SPLITLAYERS|INSANE|" << tag
                << "|col|" << ci2
                << "|origFace|" << origFace
                << "|spread|" << Foam::sqrt(d1) / max(Foam::sqrt(d0), VSMALL)
                << "|pts|" << nf << nl;
        }
    };

    forAll(columns, ci)
    {
        const Column& col = columns[ci];
        if (col.skip) { continue; }
        const label n = specs[col.si].nLayers;
        const label patchi = patches.findPatchID(specs[col.si].patchName);
        const labelList& stack = stackOf[ci];
        const labelList& assign = railsOf[railKey(col, 0)].assign;

        // rings-before-segment table
        labelList ringsBefore(col.chain.size() + 1, 0);
        for (label i = 1; i <= col.chain.size(); ++i)
        {
            label cnt = 0;
            forAll(assign, k) { if (assign[k] < i) { ++cnt; } }
            ringsBefore[i] = cnt;
        }

        // 1. wall face
        if (!faceDone.test(col.meshFacei))
        {
            meshMod.modifyFace(faces[col.meshFacei], col.meshFacei,
                               stack[0], -1, false, patchi, -1, false);
            faceDone.set(col.meshFacei);
            recordFace(col.meshFacei, "wall+rings", ci);
        }

        // 2. ring faces (reversed wall-face order -> normal away from wall)
        const face& wf = faces[col.meshFacei];
        for (label k = 0; k < n; ++k)
        {
            face rf(4);
            forAll(wf, c)
            {
                rf[c] = railsOf[railKey(col, 3 - c)].rings[k];
            }
            saneCheck(rf, col.meshFacei, ci, "ring");
            meshMod.addFace(rf, stack[k], stack[k + 1],
                            -1, -1, col.meshFacei, false, -1, -1, false);
        }

        // 3. interface faces absorbed
        for (const label fi : col.ifaces)
        {
            if (!faceDone.test(fi))
            {
                meshMod.removeFace(fi, -1);
                faceDone.set(fi);
                recordFace(fi, "ifaceRemoved", ci);
            }
        }

        // 4. top face -> remainder. Either side that is the LAST cell of a
        // surviving column maps to that column's remainder cell — this
        // also handles head-on chain meetings (top-to-top).
        if (!faceDone.test(col.topFace))
        {
            const face& tf = faces[col.topFace];
            auto remap = [&](const label celli) -> label
            {
                if (celli >= 0 && cellCol.found(celli))
                {
                    const label oc = cellCol[celli];
                    if (!columns[oc].skip
                     && columns[oc].chain.last() == celli)
                    {
                        return stackOf[oc][specs[columns[oc].si].nLayers];
                    }
                }
                return celli;
            };
            const face tfr = withForeignRings(tf);
            saneCheck(tfr, col.topFace, ci, "top");
            const label ownO = mesh.faceOwner()[col.topFace];
            if (mesh.isInternalFace(col.topFace))
            {
                const label neiO = mesh.faceNeighbour()[col.topFace];
                meshMod.modifyFace(tfr, col.topFace,
                    remap(ownO), remap(neiO), false, -1, -1, false);
            }
            else
            {
                meshMod.modifyFace(tfr, col.topFace, stack[n], -1, false,
                                   patches.whichPatch(col.topFace),
                                   -1, false);
            }
            faceDone.set(col.topFace);
            recordFace(col.topFace, "top", ci);
        }

        // 5. side faces, per chain cell, partitioned by assigned rings
        forAll(col.chain, i)
        {
            for (const label fi : cells[col.chain[i]])
            {
                if (faceDone.test(fi)) { continue; }
                if (fi == col.meshFacei || fi == col.topFace) { continue; }
                if (col.ifaces.found(fi)) { continue; }

                const face& f = faces[fi];

                // which rails does this face touch?
                FixedList<label, 4> railIdx(-1);
                label nr = 0;
                forAll(col.rails, c)
                {
                    if (f.found(col.rails[c][i])
                     && f.found(col.rails[c][i + 1]))
                    {
                        railIdx[nr++] = c;
                    }
                }
                if (nr != 2) { continue; }   // not a clean side face

                const label ca = railIdx[0], cb = railIdx[1];
                DynamicList<label> cutK;
                forAll(assign, k)
                {
                    if (assign[k] == i) { cutK.append(k); }
                }

                const label ownO = mesh.faceOwner()[fi];
                const label neiO = mesh.isInternalFace(fi)
                                 ? mesh.faceNeighbour()[fi] : -1;
                const label other =
                    (ownO == col.chain[i]) ? neiO : ownO;
                const bool otherSplit =
                    (other >= 0 && cellCol.found(other)
                  && !columns[cellCol[other]].skip);
                label otherSeg = -1;
                if (otherSplit)
                {
                    otherSeg = columns[cellCol[other]].chain.find(other);
                }

                const label nSub = cutK.size() + 1;
                for (label s = 0; s < nSub; ++s)
                {
                    auto pickPt = [&](const label c, const bool lower)
                        -> label
                    {
                        const label idx = lower ? s - 1 : s;
                        if (lower && s == 0)
                        {
                            return col.rails[c][i];
                        }
                        if (!lower && s == nSub - 1)
                        {
                            return col.rails[c][i + 1];
                        }
                        return railsOf[railKey(col, c)].rings[cutK[idx]];
                    };

                    face sf(f.size());
                    forAll(f, k)
                    {
                        label p = f[k];
                        if (p == col.rails[ca][i]) { p = pickPt(ca, true); }
                        else if (p == col.rails[ca][i+1])
                        {
                            p = pickPt(ca, false);
                        }
                        else if (p == col.rails[cb][i])
                        {
                            p = pickPt(cb, true);
                        }
                        else if (p == col.rails[cb][i+1])
                        {
                            p = pickPt(cb, false);
                        }
                        sf[k] = p;
                    }

                    const label stackIdx = ringsBefore[i] + s;
                    label sOwn = ownO, sNei = neiO;
                    if (sOwn == col.chain[i]) { sOwn = stack[stackIdx]; }
                    else if (otherSplit && sOwn == other)
                    {
                        const labelList& oStack = stackOf[cellCol[other]];
                        const labelList& oAssign =
                            railsOf[railKey(columns[cellCol[other]], 0)]
                                .assign;
                        label oBefore = 0;
                        forAll(oAssign, k)
                        {
                            if (oAssign[k] < otherSeg) { ++oBefore; }
                        }
                        sOwn = oStack[oBefore + s];
                    }
                    if (sNei == col.chain[i]) { sNei = stack[stackIdx]; }
                    else if (otherSplit && sNei == other)
                    {
                        const labelList& oStack = stackOf[cellCol[other]];
                        const labelList& oAssign =
                            railsOf[railKey(columns[cellCol[other]], 0)]
                                .assign;
                        label oBefore = 0;
                        forAll(oAssign, k)
                        {
                            if (oAssign[k] < otherSeg) { ++oBefore; }
                        }
                        sNei = oStack[oBefore + s];
                    }
                    const label sPatch = mesh.isInternalFace(fi)
                                       ? -1 : patches.whichPatch(fi);

                    const face sfr = withForeignRings(sf);
                    saneCheck(sfr, fi, ci, "strip");
                    if (s == 0)
                    {
                        meshMod.modifyFace(sfr, fi, sOwn, sNei, false,
                                           sPatch, -1, false);
                    }
                    else
                    {
                        meshMod.addFace(sfr, sOwn, sNei,
                                        -1, -1, fi, false, sPatch, -1, false);
                    }
                }
                faceDone.set(fi);
                if (debugDump)
                {
                    const label ownO2 = mesh.faceOwner()[fi];
                    const label nei2 = mesh.isInternalFace(fi)
                                     ? mesh.faceNeighbour()[fi] : -1;
                    const label oth = (ownO2 == col.chain[i]) ? nei2 : ownO2;
                    std::string t = "strips:seg" + std::to_string(i)
                        + ":cuts" + std::to_string(cutK.size());
                    if (oth < 0) { t += ":bnd"; }
                    else if (claimed[oth] >= 0)
                    {
                        const label oc = claimed[oth];
                        t += ":oc" + std::to_string(oc)
                           + (columns[oc].skip ? "SKIP" : "")
                           + ":oseg"
                           + std::to_string(columns[oc].chain.find(oth));
                    }
                    else { t += ":plain"; }
                    recordFace(fi, t.c_str(), ci);
                }
            }
        }
    }

    // ---- lateral faces of unsplit cells: insert ring points --------------
    labelHashSet lateralFaces;
    forAllConstIters(segRings, iter)
    {
        const edge& ce = iter.key();
        label meshEdgei = -1;
        for (const label ei : mesh.pointEdges()[ce[0]])
        {
            if (meshEdges[ei] == ce) { meshEdgei = ei; break; }
        }
        if (meshEdgei < 0) { continue; }
        for (const label fi : mesh.edgeFaces()[meshEdgei])
        {
            if (!faceDone.test(fi)) { lateralFaces.insert(fi); }
        }
    }

    for (const label fi : lateralFaces)
    {
        const face& f = faces[fi];
        DynamicList<label> nf(f.size() + 16);
        forAll(f, k)
        {
            nf.append(f[k]);
            const edge fe(f[k], f.nextLabel(k));
            const auto it = segRings.cfind(fe);
            if (it.good())
            {
                const labelList& ring = it.val();
                if (f[k] == it.key()[0])
                {
                    forAll(ring, j) { nf.append(ring[j]); }
                }
                else
                {
                    forAllReverse(ring, j) { nf.append(ring[j]); }
                }
            }
        }
        const label ownO = mesh.faceOwner()[fi];
        const label neiO =
            mesh.isInternalFace(fi) ? mesh.faceNeighbour()[fi] : -1;
        const face nfF{labelList(nf)};
        saneCheck(nfF, fi, -1, "lateral");
        meshMod.modifyFace
        (
            nfF, fi, ownO, neiO, false,
            mesh.isInternalFace(fi) ? -1 : patches.whichPatch(fi),
            -1, false
        );
        faceDone.set(fi);
        recordFace(fi, "lateral", -1);
    }

    // ---- apply ------------------------------------------------------------
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
    mesh.updateMesh(map());

    if (debugDump)
    {
        OFstream os(runTime.constantPath()/"splitLayersDebugMap");
        os  << "# finalFace origFace tag column" << nl;
        const labelList& fMap = map().faceMap();   // final -> origin/master
        forAll(fMap, ff)
        {
            const label of = fMap[ff];
            if (of >= 0 && faceTagOf.found(of))
            {
                os  << ff << ' ' << of << ' ' << faceTagOf[of] << ' '
                    << faceColOf[of] << nl;
            }
        }
        OFstream osc(runTime.constantPath()/"splitLayersDebugCols");
        osc << "# ci skip chain | rail0 | rail1 | rail2 | rail3" << nl;
        forAll(columns, ci)
        {
            const Column& c2 = columns[ci];
            osc << ci << ' ' << (c2.skip ? 1 : 0) << ' ';
            for (const label cc : c2.chain) { osc << cc << ' '; }
            forAll(c2.rails, r2)
            {
                osc << "| ";
                for (const label p2 : c2.rails[r2]) { osc << p2 << ' '; }
            }
            osc << nl;
        }
        Info<< "SPLITLAYERS|debugDump|written" << nl;
    }

    Info<< "SPLITLAYERS|topology|cells|" << mesh.nCells()
        << "|points|" << mesh.nPoints() << nl;

    // ---- offset-vector smoothing (wall points fixed) ----------------------
    if (nOptSweeps > 0)
    {
        pointField newPts(mesh.points());
        const labelList& pointMap = map().pointMap();
        const labelList& revPointMap = map().reversePointMap();

        Map<DynamicList<label>> children;
        forAll(pointMap, np)
        {
            const label oldp = pointMap[np];
            if (oldp >= 0 && revPointMap[oldp] != np)
            {
                children(oldp).append(np);
            }
        }

        DynamicList<label> inserted;
        DynamicList<label> insBase, insK;
        forAllConstIters(railsOf, iter)
        {
            const label b = iter.key()[0];   // rail key = (base, 1st joint)
            if (!iter.val().rings.size()) { continue; }
            if (!children.found(b)) { continue; }
            // corner-line bases carry TWO rails: their children mix both
            // rails' rings and cannot be layer-indexed by distance alone —
            // leave those unoptimized (v1)
            if (children[b].size() != iter.val().rings.size()) { continue; }
            const label newBase = revPointMap[b];
            labelList sorted(children[b]);
            std::sort(sorted.begin(), sorted.end(),
                [&](label x, label y)
                {
                    return magSqr(newPts[x] - newPts[newBase])
                         < magSqr(newPts[y] - newPts[newBase]);
                });
            forAll(sorted, j)
            {
                if (j >= iter.val().rings.size()) { break; }
                inserted.append(sorted[j]);
                insBase.append(b);
                insK.append(j);
            }
        }

        const labelListList& pointPoints = mesh.pointPoints();
        Map<label> layerOf;
        labelList baseNewOf(newPts.size(), -1);
        forAll(inserted, i)
        {
            layerOf.insert(inserted[i], insK[i]);
            baseNewOf[inserted[i]] = revPointMap[insBase[i]];
        }
        const Vector<label>& geomD = mesh.geometricD();

        for (label sweep = 0; sweep < nOptSweeps; ++sweep)
        {
            vectorField acc(newPts.size(), Zero);
            scalarField cnt(newPts.size(), 0.0);
            forAll(inserted, i)
            {
                const label np = inserted[i];
                for (const label nb : pointPoints[np])
                {
                    if (layerOf.found(nb) && layerOf[nb] == insK[i])
                    {
                        acc[np] += newPts[nb] - newPts[baseNewOf[nb]];
                        cnt[np] += 1.0;
                    }
                }
            }
            forAll(inserted, i)
            {
                const label np = inserted[i];
                if (cnt[np] < 1) { continue; }
                const point& bp = newPts[baseNewOf[np]];
                const vector vOwn = newPts[np] - bp;
                vector v = 0.5 * vOwn + 0.5 * acc[np] / cnt[np];
                const scalar mv = mag(v);
                if (mv > SMALL) { v *= mag(vOwn) / mv; }
                else { v = vOwn; }
                point target = bp + v;
                for (direction cmpt = 0; cmpt < 3; ++cmpt)
                {
                    if (geomD[cmpt] < 0) { target[cmpt] = newPts[np][cmpt]; }
                }
                newPts[np] = target;
            }
        }
        mesh.movePoints(newPts);
        Info<< "SPLITLAYERS|optimize|sweeps|" << nOptSweeps
            << "|points|" << inserted.size() << nl;
    }

    mesh.setInstance(runTime.constant());
    mesh.write();

    Info<< "SPLITLAYERS|result|split-ok" << nl << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
