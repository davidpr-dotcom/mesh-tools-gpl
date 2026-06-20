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
    label nCol = -1;            // I-INT: this column's layer count (gap-capped;
                                // == spec nLayers externally, fewer in a thin gap)
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
                       "quality-gated layer optimize sweeps — accept a smoothing move "
                       "only if neither local nonOrtho nor skew rises (default 10, 0=off)");
    argList::addOption("featureAngle", "deg",
                       "feature-point detection angle (default 45)");
    argList::addOption("convexRidgeAngle", "deg",
                       "convex-ridge geometric skip threshold, I0.1 (default 179: "
                       "the RIDGE heuristic is retired — ridges build and the "
                       "quality-retry orchestrator demotes any that fail — but "
                       "near-anti-parallel opposing walls are still skipped to mask "
                       "the opposing-wall construction bug. 180=fully off (exposes "
                       "that crash); a sharper angle re-enables the ridge skip)");
    argList::addBoolOption("dryRun", "identify and report columns only");
    argList::addBoolOption("debugDump",
        "write constant/splitLayersDebugMap: finalFace origFace tag column");
    argList::addOption("skipColumns", "file",
        "force-skip these column ids (newline list); the meshcore quality-retry "
        "orchestrator passes the columns whose built faces failed checkMesh");
    argList::addBoolOption("gapCap",
        "internal-flow gap budgeting (task 19): cap a column's layer COUNT when its "
        "full stack cannot fit the available inward depth — the oblique/curved thin "
        "gaps the anti-parallel test misses. The planner sets this ONLY for internal "
        "geometries, so external aero is byte-identical (flag off => unchanged).");
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
    const scalar convexAngle = args.getOrDefault<scalar>("convexRidgeAngle", 179.0);
    const bool dryRun = args.found("dryRun");
    const bool debugDump = args.found("debugDump");
    const bool gapCap = args.found("gapCap");   // task 19: internal-flow gap budgeting

    // Columns the meshcore quality-retry orchestrator demoted on a prior round
    // (their built faces failed checkMesh). Mesh-quality-driven, not a geometric
    // guess. Column ids are stable across runs (assigned in PASS 1, before skip).
    labelHashSet forcedSkip;
    {
        fileName scf;
        if (args.readIfPresent("skipColumns", scf))
        {
            std::ifstream scs(scf.c_str());
            label c;
            while (scs >> c) { forcedSkip.insert(c); }
            Info<< "SPLITLAYERS|skipColumns|count|" << forcedSkip.size() << nl;
        }
    }

    // forensics: original-face -> (construct tag, column id)
    Map<word> faceTagOf;
    Map<label> faceColOf;
    // task 19 (opposing-wall gap collision): for a strip side face, the OTHER
    // column claiming the cell across the face (claimed[oth]). Lets the OPENCELL
    // attributor name BOTH facing-wall columns and dump their seg0 ring geometry.
    Map<label> faceOtherColOf;
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
    // I-INT G1 (emit-only, NO behaviour change): the gap budget. How many spec
    // layers fit the available rail depth at FULL size (cumHeights[k] <=
    // minRailLen)? The rail stops where it runs out of room — the opposing wall,
    // or the gap centre where the opposing column's claim blocks it — so
    // minRailLen IS the available depth: ~half-gap internally, ~needed externally.
    // External columns fit all layers; thin-gap columns fit fewer. This proves the
    // internal/external discriminator before any topology change (G2).
    const scalar gapMargin = 0.7;   // thin-gap: layers fill <= 70% of avail depth,
                                    // leaving a healthy remainder + meeting buffer
    List<label> gapCapped(specs.size(), 0);          // cols with nCol < nLayers
    List<label> gapSkip(specs.size(), 0);            // cols skipped (gap < 1 layer)
    List<label> gapTot(specs.size(), 0);             // total built cols
    List<label> gapMinFit(specs.size(), labelMax);   // min nCol over cols
    // task 19 confirmatory (emit-only): count + sample columns that keep FULL
    // layers while crammed into a thin space (scale<0.5) — the open-cell danger
    // zone — and dump their gap-discriminator decision.
    label gapDanger = 0;

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

            label blockCol = -1;    // I-INT G3: column whose claim stopped the rail
            label blockFace = -1;   // or the opposing-patch wall face it reached
            while (minRailLen() < needed[si] && mesh.isInternalFace(iface))
            {
                const label cOwn = mesh.faceOwner()[iface];
                const label cNei = mesh.faceNeighbour()[iface];
                const label next = (cOwn == col.chain.last()) ? cNei : cOwn;
                if (next < 0) { break; }
                if (claimed[next] >= 0) { blockCol = claimed[next]; break; }
                if (patchFaceCount[next] > 0)
                {
                    // rail reached a boundary cell — capture its face on THIS
                    // layered patch (the opposing wall of a thin gap).
                    forAll(cells[next], fj)
                    {
                        const label bf = cells[next][fj];
                        if (!mesh.isInternalFace(bf)
                         && patches.whichPatch(bf) == pp.index())
                        { blockFace = bf; break; }
                    }
                    break;
                }
                if (!isHex(next)) { break; }
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
            // I-INT G2: gap-aware layer count. minRailLen is the available depth.
            // If ALL spec layers fit at full size -> external, keep the exact old
            // scale (byte-identical). Else (thin gap) -> as many FULL-size layers
            // as fit within gapMargin·avail (full y1 preserved, healthy remainder);
            // 0 -> skip (gap too thin for even one layer, cleaner than degenerate).
            {
                const scalar avail = minRailLen();
                // I-INT G3: cap ONLY at a TRUE gap — an ANTI-PARALLEL opposing
                // wall (rail blocked by the opposing column's claim, or reaching
                // that wall directly). A reentrant CORNER's blocker is ~perpendic-
                // ular (dot ~0) so it is NOT capped -> corners/externals stay
                // byte-identical. -0.5 = walls facing within 60° of head-on.
                vector wn = mesh.faceAreas()[col.meshFacei];
                wn /= (mag(wn) + VSMALL);
                label oppFace = -1;
                if (blockCol >= 0 && blockCol < columns.size())
                {
                    oppFace = columns[blockCol].meshFacei;   // opposing col's wall
                }
                else if (blockFace >= 0)
                {
                    oppFace = blockFace;                     // opposing wall reached
                }
                else if (!mesh.isInternalFace(col.topFace)
                      && patches.whichPatch(col.topFace) == pp.index())
                {
                    oppFace = col.topFace;
                }
                bool isGap = false;
                if (oppFace >= 0)
                {
                    vector on = mesh.faceAreas()[oppFace];
                    on /= (mag(on) + VSMALL);
                    isGap = ((wn & on) < -0.5);
                }

                // task 19 fix: oblique/curved thin gaps. The anti-parallel test
                // (wn.on < -0.5) misses gaps whose opposing wall is oblique — the
                // manifold's are wn.on +0.1..+0.4 (GAPDECIDE, 2026-06-17), so 98.5%
                // of its columns kept full layers crammed into a sub-mm depth. For
                // an INTERNAL-flow geometry the planner sets -gapCap; then ANY
                // column whose full stack cannot fit the available inward depth is
                // gap-constrained, so it takes the cap branch below (nFit full-size
                // layers within avail*gapMargin, leaving a healthy remainder) — no
                // dependence on the opposing-wall angle. External aero never sets
                // the flag, so it is byte-identical (this block is a no-op there).
                if (gapCap && cumHeights[si].last() > avail + SMALL)
                {
                    isGap = true;
                }

                if (!isGap || cumHeights[si].last() <= avail + SMALL)
                {
                    col.nCol = specs[si].nLayers;
                    col.scale = min(scalar(1), avail / needed[si]);
                }
                else
                {
                    label nFit = 0;
                    forAll(cumHeights[si], k)
                    {
                        if (cumHeights[si][k] <= avail*gapMargin) { ++nFit; }
                    }
                    col.nCol = nFit;
                    col.scale = scalar(1);   // full y1; the nFit layers fit by build
                    if (nFit == 0)
                    {
                        col.skip = true; col.reason = "thinGap"; ++gapSkip[si];
                    }
                }
                ++gapTot[si];
                if (col.nCol < specs[si].nLayers) { ++gapCapped[si]; }
                gapMinFit[si] = min(gapMinFit[si], col.nCol);

                // task 19 confirmatory (emit-only): a column with FULL layers
                // crammed into a thin space (scale<0.5) is the open-cell danger
                // zone. Dump the gap-discriminator decision to see WHY isGap
                // wasn't set (oppFace missing? wn.on not anti-parallel?) and to
                // confirm the remainder is ~0 (avail ~= scale*cumLast).
                if (col.nCol == specs[si].nLayers && col.scale < 0.5)
                {
                    ++gapDanger;
                    if (gapDanger <= 60)
                    {
                        scalar won = 2.0;     // sentinel: no oppFace found
                        label oppPatch = -1;
                        if (oppFace >= 0)
                        {
                            vector onv = mesh.faceAreas()[oppFace];
                            onv /= (mag(onv) + VSMALL);
                            won = (wn & onv);
                            oppPatch = mesh.isInternalFace(oppFace)
                                     ? -1 : patches.whichPatch(oppFace);
                        }
                        Info<< "SPLITLAYERS|GAPDECIDE|col|" << columns.size()
                            << "|chain|" << col.chain.size()
                            << "|blockCol|" << blockCol
                            << "|blockFace|" << blockFace
                            << "|topInternal|" << mesh.isInternalFace(col.topFace)
                            << "|oppFace|" << oppFace
                            << "|oppPatch|" << oppPatch
                            << "|wn.on|" << won
                            << "|avail|" << avail
                            << "|needed|" << needed[si]
                            << "|cumLast|" << cumHeights[si].last()
                            << "|isGap|" << isGap
                            << "|nCol|" << col.nCol
                            << "|scale|" << col.scale << nl;
                    }
                }
            }
            columns.append(col);
        }
    }

    forAll(specs, si)
    {
        Info<< "SPLITLAYERS|gap|" << specs[si].patchName
            << "|cappedCols|" << gapCapped[si]
            << "|skippedThinGap|" << gapSkip[si]
            << "|totalCols|" << gapTot[si]
            << "|minFit|" << (gapTot[si] ? gapMinFit[si] : specs[si].nLayers)
            << "|nLayers|" << specs[si].nLayers << nl;
    }
    Info<< "SPLITLAYERS|gapDangerTotal|" << gapDanger << nl;

    // ---- task 19: gradual layer-count termination (gapCap only) ----------
    // The count cap leaves abrupt nCol jumps (e.g. 8 next to 2) where a column
    // that fit its full stack abuts a pinched neighbour. The variable-n strip/
    // lateral topology only closes GENTLE steps (the `term` clamp handles +-1),
    // so big jumps dangle (BADFACE|lateral, ownN/neiN = 8/2, 2026-06-17). Smooth
    // the count field so edge-adjacent SPLIT columns differ by <= maxNStep: an
    // 8->2 cliff becomes an 8->..->2 gradient every transition can close.
    // Monotone-decreasing => terminates; the fixpoint is order-independent
    // (min-plus relaxation). Split-vs-skip (nCol 0) transitions are left alone
    // (the lateral pass handles them; they do not dangle). External aero never
    // sets -gapCap => this whole block is a no-op there (byte-identical).
    if (gapCap)
    {
        const label maxNStep = 1;
        label nLowered = 0, totSweeps = 0;
        forAll(specs, si)
        {
            const polyPatch& pp =
                patches[patches.findPatchID(specs[si].patchName)];
            labelList colOf(pp.size(), -1);             // patch face -> column id
            forAll(columns, ci)
            {
                const Column& c = columns[ci];
                if (c.si == si && !c.skip && c.nCol > 0)
                {
                    colOf[c.patchFacei] = ci;
                }
            }
            const labelListList& ff = pp.faceFaces();   // edge-connected patch faces
            bool changed = true;
            label sweeps = 0;
            while (changed && sweeps < 1000)
            {
                changed = false;
                ++sweeps;
                forAll(colOf, pfi)
                {
                    const label ci = colOf[pfi];
                    if (ci < 0) { continue; }
                    label lo = labelMax;
                    forAll(ff[pfi], k)
                    {
                        const label nci = colOf[ff[pfi][k]];
                        if (nci >= 0) { lo = min(lo, columns[nci].nCol); }
                    }
                    if (lo != labelMax && columns[ci].nCol > lo + maxNStep)
                    {
                        columns[ci].nCol = lo + maxNStep;
                        ++nLowered;
                        changed = true;
                    }
                }
            }
            totSweeps += sweeps;
        }
        Info<< "SPLITLAYERS|gapSmooth|maxStep|" << maxNStep
            << "|lowered|" << nLowered << "|sweeps|" << totSweeps << nl;
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

    // ---- S1: min-scale consensus (I1.1) ---------------------------------
    // Columns that share a GENUINE rail (same first edge AND identical full
    // joint path) but carry different per-column scales — because their OTHER
    // rails differ in length — trip the scale arm of the conflict check, so the
    // later column is skipped (the "scaleOnly" conflicts; S0: bracket 31,
    // manifold 221k). Union such columns and give the component ONE scale (its
    // MIN, so no column claims more height than it can reach); the existing
    // first-come conflict check then passes (equal scale -> equal assign) and
    // the layers are kept. Strictly conflict-reducing (joints/nSeg untouched).
    //
    // NOTE (2026-06-16): broadening this union to share-first-edge prefix-sharers
    // (to lift S2b recovery) was tried and REVERTED — it cascades. A column needs
    // ONE scale for its 4 rails' ring faces to stay planar, and rails are shared
    // broadly (even flat structured regions), so harmonizing scale across shared
    // rails couples the whole connected surface into one component (measured:
    // bracket nComp 1, meanScale 0.146 — the BL at 14.6% everywhere). Lifting the
    // remaining corner bisectors is the I1.3 corner-column problem, not a scale
    // tweak. The consensus|nComp diagnostic below monitors for any such cascade.
    {
        labelList parent(columns.size());
        forAll(parent, i) { parent[i] = i; }
        auto root = [&](label a) -> label   // iterative find + path-halving
        {
            while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
            return a;
        };
        auto unite = [&](label a, label b)
        {
            const label ra = root(a), rb = root(b);
            if (ra != rb) { parent[max(ra, rb)] = min(ra, rb); }
        };
        // group rail instances by first-edge key; within a bucket, columns whose
        // rail has the IDENTICAL joint path share a genuine rail -> union them.
        EdgeMap<DynamicList<label>> byKey;     // key -> packed (ci*4 + railIndex)
        forAll(columns, ci)
        {
            if (columns[ci].skip || forcedSkip.found(ci)) { continue; }
            forAll(columns[ci].rails, c)
            {
                DynamicList<label>& bucket = byKey(railKey(columns[ci], c));
                for (const label packed : bucket)
                {
                    if (columns[packed/4].rails[packed%4] == columns[ci].rails[c])
                    {
                        unite(ci, packed/4);
                    }
                }
                bucket.append(ci*4 + c);
            }
        }
        Map<scalar> compMin;
        forAll(columns, ci)
        {
            if (columns[ci].skip || forcedSkip.found(ci)) { continue; }
            const label r = root(ci);
            if (!compMin.found(r) || columns[ci].scale < compMin[r])
            {
                compMin.set(r, columns[ci].scale);
            }
        }
        label nLowered = 0;
        forAll(columns, ci)
        {
            if (columns[ci].skip || forcedSkip.found(ci)) { continue; }
            const scalar s = compMin[root(ci)];
            if (s < columns[ci].scale - SMALL) { ++nLowered; }
            columns[ci].scale = s;
        }
        Info<< "SPLITLAYERS|consensus|scaleLowered|" << nLowered << nl;
        // diagnostic (behaviour-neutral): is the union cascading, and how hard is
        // the resulting compression? maxComp == ~nCols => one giant component (a
        // cascade); minScale/meanScale say whether that actually thins the BL.
        {
            Map<label> compSize;
            scalar minScale = GREAT, sumScale = 0;
            label nCols = 0;
            forAll(columns, ci)
            {
                if (columns[ci].skip || forcedSkip.found(ci)) { continue; }
                ++compSize(root(ci));
                minScale = min(minScale, columns[ci].scale);
                sumScale += columns[ci].scale;
                ++nCols;
            }
            label maxComp = 0;
            forAllConstIters(compSize, it) { maxComp = max(maxComp, it.val()); }
            Info<< "SPLITLAYERS|consensus|nComp|" << compSize.size()
                << "|maxComp|" << maxComp
                << "|minScale|" << minScale
                << "|meanScale|" << (nCols ? sumScale/nCols : 1.0) << nl;
        }
    }

    List<label> nConflict(specs.size(), 0), nWarp(specs.size(), 0),
                nCross(specs.size(), 0);
    List<label> nConvex(specs.size(), 0);        // skippedConvexRidge (I0.1)
    // S2 diagnostic (I1.1, behaviour-neutral): for joints-arm conflicts (corner
    // bisectors), histogram the shared-prefix depth (how many leading rail points
    // the two diverging rails share). Depth ~2 = share only seg0 → per-segment
    // consensus recovers; deeper = true corner → I1.3. Decides the S2b approach.
    Map<label> jointsPrefixHist;

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
        // I-INT: only this column's (gap-capped) layers, not the full spec count
        const label nc = (col.nCol >= 0) ? col.nCol : cum.size();
        labelList assign(nc);
        for (label k = 0; k < nc; ++k)
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

        // forced skip from the meshcore quality-retry orchestrator: this column's
        // built faces failed checkMesh on a prior round (mesh-quality-driven, not
        // a geometric guess). Demote it; its neighbours may recover.
        if (forcedSkip.found(ci))
        {
            col.skip = true; col.reason = "forcedSkip";
            continue;
        }

        // --- convex-ridge opposed-column skip (I0.1 cow ears, 2026-06-13) ---
        // Two adjacent split columns sharing a sharp CONVEX wall edge get rails
        // that diverge over the bend; their seg0 strips meet as strips:opposed
        // and come out wrong-oriented/highly-skew. As a RIDGE heuristic this is
        // RETIRED: the default convexRidgeAngle is 179, so genuine ridges (cow
        // ears, bracket corners) are NO LONGER skipped — they build and the
        // meshcore quality-retry orchestrator demotes only the columns whose faces
        // actually fail checkMesh (mesh-quality-driven, not a geometric guess).
        // At 179 this block STILL fires for near-anti-parallel opposing walls
        // (e.g. a narrow gap): those columns hit a separate opposing-wall
        // construction bug (polyTopoChange::getFaceOrder, boundary face w/o patch)
        // and must be skipped until that bug is fixed (then default -> 180 = fully
        // off). 180 disables this block entirely; a sharper angle (e.g. 70) brings
        // back the full ridge heuristic. Real corner fix = I1.3/I1.4 (warp cols).
        if (convexAngle < 180.0)
        {
            const label patchi   = patches.findPatchID(specs[col.si].patchName);
            const label wf       = col.meshFacei;                 // wall face
            const vector na      = mesh.faceAreas()[wf];
            const vector nahat   = na/(mag(na) + VSMALL);          // outward of domain
            const point  Ca      = mesh.faceCentres()[wf];
            const scalar cosConvex = Foam::cos(degToRad(convexAngle));
            bool convexRidge = false;
            for (const label ei : mesh.faceEdges()[wf])
            {
                for (const label nb : mesh.edgeFaces()[ei])
                {
                    if (nb == wf || mesh.isInternalFace(nb))   continue;
                    if (patches.whichPatch(nb) != patchi)      continue;
                    const vector nbA   = mesh.faceAreas()[nb];
                    const vector nbhat = nbA/(mag(nbA) + VSMALL);
                    if ((nahat & nbhat) >= cosConvex)         continue;  // not sharp enough
                    const point  Cb    = mesh.faceCentres()[nb];
                    if ((nahat & (Cb - Ca)) <= 0)             continue;  // concave -> not ours
                    const label nbCell = mesh.faceOwner()[nb];
                    if (claimed[nbCell] >= 0 && claimed[nbCell] != ci)
                    {
                        convexRidge = true; break;            // adjacent claimed column on a convex ridge
                    }
                }
                if (convexRidge) break;
            }
            if (convexRidge)
            {
                col.skip = true; col.reason = "convexRidge";
                ++nConvex[col.si];
                continue;
            }
        }

        const label m = col.chain.size();

        // (a) shared-rail consistency: chain length, scale AND the ring->
        // segment assignment table must match what a neighbouring column
        // already recorded on a shared rail. Assignment mismatch was the
        // bracket quality defect (2026-06-12 forensics: all 89 bad faces
        // were strips of columns whose shared rails were subdivided under
        // a neighbour's table while their own rails used their own).
        const labelList a0 = computeAssign(col, 0);
        bool conflict = false;
        label jointsPrefix = -1;   // S2 diag: shallowest shared-prefix on a joints-mismatch rail
        forAll(col.rails, c)
        {
            const auto it = railsOf.cfind(railKey(col, c));
            if (!it.good()) { continue; }
            const bool jb = (it.val().joints != col.rails[c]);   // FULL-PATH identity
            if (it.val().nSeg != m
             || mag(it.val().scale - col.scale) > 1e-6 * col.scale + VSMALL
             || it.val().assign != a0
             || jb)
            {
                conflict = true;
            }
            if (jb)   // common leading points before the two rails diverge
            {
                const labelList& A = it.val().joints;
                const labelList& B = col.rails[c];
                const label n = min(A.size(), B.size());
                label p = 0;
                while (p < n && A[p] == B[p]) { ++p; }
                if (jointsPrefix < 0 || p < jointsPrefix) { jointsPrefix = p; }
            }
        }
        if (conflict)
        {
            col.skip = true; col.reason = "conflict";
            ++nConflict[col.si];
            if (jointsPrefix >= 0) { jointsPrefixHist(jointsPrefix) += 1; }  // S2 diag
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
                << "|skippedConvexRidge|" << nConvex[si]
                << nl;
        }
        // S2 diagnostic: shared-prefix-depth histogram of the joints conflicts.
        Info<< "SPLITLAYERS|jointsConflict|prefixHistogram";
        forAllConstIters(jointsPrefixHist, it)
        {
            Info<< "|" << it.key() << ":" << it.val();
        }
        Info<< nl;
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

    // ---- I1.3 (I3.1): convex-ridge wall points --------------------------
    // At a sharp CONVEX wall edge, adjacent columns grow DIVERGENT rails from its
    // endpoints (different first-interior joints -> different railKeys -> separate
    // ring sets), so their seg0 strips meet OPPOSED and fail to close (the I0.1
    // cow-ear / manifold open-cell defect). Mark those endpoints; I3.2 gives them
    // ONE consensus ring column so the strips coincide. Concave corners are
    // excluded by the face-centre test (same as the retired convex skip).
    // (Emit-only here — proves the detector before any construction change.)
    labelHashSet ridgePts;
    {
        const scalar cosRidge = Foam::cos(degToRad(featAngle));
        forAll(specs, si)
        {
            const label patchi = patches.findPatchID(specs[si].patchName);
            if (patchi < 0) { continue; }
            const polyPatch& pp = patches[patchi];
            forAll(pp, pfi)
            {
                const label wf = pp.start() + pfi;
                const vector na = mesh.faceAreas()[wf];
                const vector nahat = na/(mag(na) + VSMALL);
                const point Ca = mesh.faceCentres()[wf];
                for (const label ei : mesh.faceEdges()[wf])
                {
                    for (const label nb : mesh.edgeFaces()[ei])
                    {
                        if (nb == wf || mesh.isInternalFace(nb)) { continue; }
                        if (patches.whichPatch(nb) != patchi) { continue; }
                        const vector nbA = mesh.faceAreas()[nb];
                        const vector nbhat = nbA/(mag(nbA) + VSMALL);
                        if ((nahat & nbhat) >= cosRidge) { continue; } // not sharp
                        const point Cb = mesh.faceCentres()[nb];
                        if ((nahat & (Cb - Ca)) <= 0) { continue; }    // concave
                        const edge& e = mesh.edges()[ei];
                        ridgePts.insert(e[0]);
                        ridgePts.insert(e[1]);
                    }
                }
            }
        }
        Info<< "SPLITLAYERS|ridgePts|" << ridgePts.size() << nl;
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

            rd.rings.setSize(col.nCol);   // I-INT: this column's gap-capped count
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

    // S2b — recovery store for prefix-divergent (corner-bisector) sharers,
    // populated by the recovery pass. EMPTY here, so the accessors below fall
    // through to the canonical railsOf and behaviour is byte-identical to today;
    // the recovery increment fills these with a sharer's OWN rail rings (reused
    // canonical prefix + own divergent suffix) and assign table.
    Map<FixedList<labelList, 4>> sharerRings;   // ci -> 4 rails' recovered rings
    Map<labelList> sharerAssign;                // ci -> recovered assign

    // ---- S2b: recover prefix-divergent (corner-bisector) sharers ----------
    // A column skipped as "conflict" because a rail shares a canonical's first
    // edge then DIVERGES (S2 diagnostic: every bracket/cow joints conflict is a
    // shared prefix). Rebuild it: adopt the (min) canonical scale; on each rail
    // REUSE the canonical's shared-prefix ring LABELS (so the shared edges stay
    // single-valued) and create OWN points on the divergent suffix. Every step
    // is guarded — scale must match to reuse, the column must be warp-free, and
    // a suffix segment claimed by another base aborts the recovery — so anything
    // that cannot recover cleanly STAYS SKIPPED (coverage-only fallback). The
    // canonical creation path above is UNTOUCHED.
    {
        label nRecovered = 0;
        // S2b abort-reason census (behaviour-neutral): how many conflict columns
        // abort, and on which guard — sizes whether a scale-harmonization lift is
        // worth building vs handing the rest to I1.3 corner columns.
        label nNoCanon = 0, nWarp = 0, nScale = 0, nEdgeClaim = 0;
        forAll(columns, ci)
        {
            Column& col = columns[ci];
            if (!col.skip || col.reason == nullptr) { continue; }
            if (std::string(col.reason) != "conflict") { continue; }

            // adopt the MIN canonical scale over this column's colliding rails
            scalar s = col.scale;
            bool anyCanon = false;
            forAll(col.rails, c)
            {
                const auto it = railsOf.cfind(railKey(col, c));
                if (it.good() && it.val().rings.size())
                {
                    s = min(s, it.val().scale);
                    anyCanon = true;
                }
            }
            if (!anyCanon) { ++nNoCanon; continue; }
            const scalar saved = col.scale;
            col.scale = s;                          // assign/heights now consistent at s

            // warp guard: all four rails must agree on ring->segment assignment
            const labelList a0 = computeAssign(col, 0);
            bool ok = true;
            const char* why = nullptr;
            for (label c = 1; c < 4 && ok; ++c)
            {
                if (computeAssign(col, c) != a0) { ok = false; why = "warp"; }
            }

            const label nL = col.nCol;   // I-INT: gap-capped layer count
            FixedList<labelList, 4> myRings;        // final per-layer ring labels
            FixedList<List<point>, 4> createPos;    // planned positions for to-create slots

            // PLAN (no mesh side effects): source each rail's rings PER SEGMENT
            // from segRings — the SAME registry the lateral faces read. Reuse a
            // segment's registered labels if present (so ring faces and lateral
            // faces can never disagree on a shared edge); else plan a creation.
            // Any guard failure aborts the WHOLE recovery before a single point is
            // added, so an aborted recovery leaves zero orphaned points.
            forAll(col.rails, c)
            {
                if (!ok) { break; }
                const labelList& r = col.rails[c];
                myRings[c].setSize(nL, -1);
                createPos[c].setSize(nL, point(Zero));
                scalarList arc(r.size(), Zero);
                for (label j = 1; j < r.size(); ++j)
                {
                    arc[j] = arc[j-1] + mag(pts[r[j]] - pts[r[j-1]]);
                }
                for (label j = 0; j + 1 < r.size() && ok; ++j)
                {
                    // layers (in k order) that fall on this segment
                    DynamicList<label> ks;
                    forAll(a0, k) { if (a0[k] == j) { ks.append(k); } }
                    if (ks.empty()) { continue; }
                    const edge se(r[j], r[j+1]);
                    const auto sr = segRings.cfind(se);
                    if (sr.good())
                    {
                        // REUSE the registered labels — only if the layer set
                        // matches (same scale+arc => same count); a mismatch is
                        // irreconcilable -> fall back to the strict skip.
                        if (sr.val().size() != ks.size())
                        { ok = false; why = "scale"; break; }
                        forAll(ks, i) { myRings[c][ks[i]] = sr.val()[i]; }
                    }
                    else
                    {
                        // CREATE: the edge must be free or already this base's.
                        const auto ec = edgeClaim.cfind(se);
                        if (ec.good() && ec.val() != r[0])
                        { ok = false; why = "edgeClaim"; break; }
                        forAll(ks, i)
                        {
                            const label k = ks[i];
                            const scalar h = s * cumHeights[col.si][k];
                            const scalar f =
                                (h - arc[j]) / max(arc[j+1] - arc[j], VSMALL);
                            createPos[c][k] =
                                pts[r[j]] + f*(pts[r[j+1]] - pts[r[j]]);
                            // myRings[c][k] stays -1 -> filled at commit
                        }
                    }
                }
            }
            if (!ok)
            {
                if (why == nullptr || std::string(why) == "scale") { ++nScale; }
                else if (std::string(why) == "warp") { ++nWarp; }
                else { ++nEdgeClaim; }
                col.scale = saved;                  // revert; stay skipped (fallback)
                continue;
            }

            // COMMIT (the only mesh side effects, now every guard has passed):
            // create the planned points, then register + claim each CREATED
            // segment so segRings (lateral faces) and myRings (ring faces) carry
            // identical labels on every edge of this column.
            forAll(col.rails, c)
            {
                const labelList& r = col.rails[c];
                forAll(myRings[c], k)
                {
                    if (myRings[c][k] == -1)
                    {
                        const point& rp = createPos[c][k];
                        myRings[c][k] = meshMod.addPoint(rp, r[0], -1, true);
                        newPtPos.insert(myRings[c][k], rp);
                    }
                }
                for (label j = 0; j + 1 < r.size(); ++j)
                {
                    const edge se(r[j], r[j+1]);
                    if (segRings.found(se)) { continue; }   // reused — keep labels
                    DynamicList<label> inSeg;
                    forAll(a0, k) { if (a0[k] == j) { inSeg.append(myRings[c][k]); } }
                    if (inSeg.size())
                    {
                        segRings.insert(se, labelList(inSeg));
                        edgeClaim.insert(se, -1);
                    }
                }
            }
            sharerRings.insert(ci, myRings);
            sharerAssign.insert(ci, a0);
            col.skip = false;
            col.reason = nullptr;
            ++nRecovered;
        }
        Info<< "SPLITLAYERS|recovered|" << nRecovered << nl;
        Info<< "SPLITLAYERS|recoverAbort|noCanon|" << nNoCanon
            << "|warp|" << nWarp << "|scale|" << nScale
            << "|edgeClaim|" << nEdgeClaim << nl;
    }

    // ---- PASS 4: topology -------------------------------------------------
    Map<labelList> stackOf;     // column id -> n+1 new cell labels
    Map<label> cellCol;         // chain cell -> column id
    forAll(columns, ci)
    {
        const Column& col = columns[ci];
        if (col.skip) { continue; }
        const label n = col.nCol;   // I-INT: gap-capped layer count
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

    // S2b ring/assign access — a column's OWN rail data: a recovered sharer's
    // reused-prefix+own-suffix rings (sharerRings/sharerAssign), else the
    // canonical railsOf. With the store empty this is exactly railsOf (identity).
    auto ringsFor = [&](const label ci2, const Column& cl, const label c)
        -> const labelList&
    {
        const auto sit = sharerRings.cfind(ci2);
        return sit.good() ? sit.val()[c] : railsOf[railKey(cl, c)].rings;
    };
    auto assignFor = [&](const label ci2, const Column& cl) -> const labelList&
    {
        const auto sit = sharerAssign.cfind(ci2);
        return sit.good() ? sit.val() : railsOf[railKey(cl, 0)].assign;
    };

    // OOB instrument (find the UB behind the changeMesh boundary-face crash): the
    // lateral pass indexes stack[ringsBefore[i]+s] and oStack[oBefore+s]; for a
    // recovered column whose assign doesn't line up with a neighbour these can run
    // past the (nLayers+1)-entry stack — an out-of-bounds List read that returns a
    // garbage neighbour label and crashes changeMesh. Detect + report the
    // site/column/index, and CLAMP so the run completes (mesh wrong but no UB).
    label nOOB = 0;
    auto chk = [&](const label idx, const label sz,
                   const char* site, const label cc) -> label
    {
        if (idx >= 0 && idx < sz) { return idx; }
        ++nOOB;
        if (nOOB <= 12)
        {
            Info<< "SPLITLAYERS|OOB|" << site << "|ci|" << cc
                << "|idx|" << idx << "|size|" << sz << nl;
        }
        return (idx < 0) ? 0 : sz - 1;
    };
    // I-INT layer termination: when THIS column has more layers than a NEIGHBOUR
    // (variable-n gap caps), its upper layers connect to the neighbour's REMAINDER
    // (last stack cell). This is by design, so clamp SILENTLY — unlike chk, which
    // alarms because an over-run of a column's OWN stack would be a real bug.
    auto term = [&](const label idx, const label sz) -> label
    {
        return (idx >= sz) ? sz - 1 : ((idx < 0) ? 0 : idx);
    };

    // BADFACE instrument: the changeMesh crash is a dangling face that ends as
    // (owner, neighbour -1, region -1) — a boundary face with no patch. Catch the
    // exact condition at the source (a face being modified/added with nei<0 AND
    // patch<0) and report which pass/column/face + the owner-resolution context.
    label nBad = 0, nDefer = 0;
    // cellCol holds every REMOVED chain cell (filled in the PASS-4 stack loop
    // below). A face whose owner/neighbour resolves to a cell still in cellCol
    // dangles after changeMesh (that cell is gone) — so does (nei<0 && patch<0).
    auto badFace = [&](const char* pass, const label fi, const label cc,
                       const label own, const label nei, const label patch,
                       const label other2, const bool otherSplit2)
    {
        const bool danglOwn = (own >= 0 && cellCol.found(own));
        const bool danglNei = (nei >= 0 && cellCol.found(nei));
        if ((nei < 0 && patch < 0) || danglOwn || danglNei)
        {
            ++nBad;
            if (nBad <= 16)
            {
                Info<< "SPLITLAYERS|BADFACE|" << pass << "|fi|" << fi
                    << "|ci|" << cc << "|own|" << own << "|nei|" << nei
                    << "|patch|" << patch << "|other|" << other2
                    << "|otherSplit|" << otherSplit2
                    << "|danglOwn|" << danglOwn << "|danglNei|" << danglNei
                    // task 19: name the columns owning the dangling sides + their
                    // nCol, to read the variable-n jump across the failing face
                    // (asymmetric ownN!=neiN => a deeper-neighbour defer not picked
                    // up; equal => a different lateral-handling gap).
                    << "|ownCol|" << (danglOwn ? cellCol[own] : -1)
                    << "|ownN|" << (danglOwn ? columns[cellCol[own]].nCol : -1)
                    << "|neiCol|" << (danglNei ? cellCol[nei] : -1)
                    << "|neiN|" << (danglNei ? columns[cellCol[nei]].nCol : -1)
                    << "|internalFi|" << mesh.isInternalFace(fi)
                    << "|whichPatch|" << patches.whichPatch(fi) << nl;
            }
        }
    };

    forAll(columns, ci)
    {
        const Column& col = columns[ci];
        if (col.skip) { continue; }
        const label n = col.nCol;   // I-INT: gap-capped layer count
        const label patchi = patches.findPatchID(specs[col.si].patchName);
        const labelList& stack = stackOf[ci];
        const labelList& assign = assignFor(ci, col);

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
                rf[c] = ringsFor(ci, col, 3 - c)[k];
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
                        return stackOf[oc][columns[oc].nCol];   // I-INT: remainder
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
                // If a side of this top face is the NON-last chain cell of ANOTHER
                // splitting column, the face is that column's SIDE face — its
                // correct stack-cell assignment belongs to that column's side-face
                // pass (section 5 below), which reassigns BOTH sides via the strip
                // logic. remap here can only reach the remainder, which is wrong
                // for an offset opposing-wall gap meeting (open cells) AND leaves
                // the removed cell dangling (the crash). DEFER: leave it unhandled
                // so the owning side-face pass picks it up. Head-on last-to-last
                // is NOT a side-of-other, so it still maps via remap -> remainder.
                auto sideOfOther = [&](const label c) -> bool
                {
                    if (c < 0 || !cellCol.found(c)) { return false; }
                    const label oc = cellCol[c];
                    return oc != ci && !columns[oc].skip
                        && columns[oc].chain.last() != c;
                };
                if (sideOfOther(ownO) || sideOfOther(neiO))
                {
                    ++nDefer;       // owning column's side-face pass handles it
                }
                else
                {
                    meshMod.modifyFace(tfr, col.topFace,
                        remap(ownO), remap(neiO), false, -1, -1, false);
                    faceDone.set(col.topFace);
                    recordFace(col.topFace, "top", ci);
                }
            }
            else
            {
                meshMod.modifyFace(tfr, col.topFace, stack[n], -1, false,
                                   patches.whichPatch(col.topFace),
                                   -1, false);
                faceDone.set(col.topFace);
                recordFace(col.topFace, "top", ci);
            }
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

                // I-INT G3 Step 2: variable-n layer termination. If the neighbour
                // is DEEPER (more layers), it must build this shared face with ITS
                // full cuts so all its layers get side faces; our shallower side
                // terminates to our remainder via the `term` clamp when it does.
                // DEFER — leave it for the deeper column's side-face pass (which
                // also iterates this face). Equal counts -> first-come, consistent;
                // shallower-neighbour -> WE build and IT terminates against us.
                if (otherSplit && columns[cellCol[other]].nCol > col.nCol)
                {
                    ++nDefer;
                    continue;
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
                        return ringsFor(ci, col, c)[cutK[idx]];
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
                    if (sOwn == col.chain[i])
                    { sOwn = stack[chk(stackIdx, stack.size(), "stack:sOwn", ci)]; }
                    else if (otherSplit && sOwn == other)
                    {
                        const labelList& oStack = stackOf[cellCol[other]];
                        const labelList& oAssign =
                            assignFor(cellCol[other], columns[cellCol[other]]);
                        label oBefore = 0;
                        forAll(oAssign, k)
                        {
                            if (oAssign[k] < otherSeg) { ++oBefore; }
                        }
                        sOwn = oStack[term(oBefore + s, oStack.size())];
                    }
                    if (sNei == col.chain[i])
                    { sNei = stack[chk(stackIdx, stack.size(), "stack:sNei", ci)]; }
                    else if (otherSplit && sNei == other)
                    {
                        const labelList& oStack = stackOf[cellCol[other]];
                        const labelList& oAssign =
                            assignFor(cellCol[other], columns[cellCol[other]]);
                        label oBefore = 0;
                        forAll(oAssign, k)
                        {
                            if (oAssign[k] < otherSeg) { ++oBefore; }
                        }
                        sNei = oStack[term(oBefore + s, oStack.size())];
                    }
                    const label sPatch = mesh.isInternalFace(fi)
                                       ? -1 : patches.whichPatch(fi);

                    const face sfr = withForeignRings(sf);
                    saneCheck(sfr, fi, ci, "strip");
                    badFace("strip", fi, ci, sOwn, sNei, sPatch,
                            other, otherSplit);
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
                        faceOtherColOf.insert(fi, oc);   // task 19: opposing column
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

    // task 19 never-crash guard: a lateral face between two SPLIT columns can
    // escape the side-face pass (a transition defer not picked up — ~40 residual
    // after gapSmooth) and arrive here with a REMOVED chain cell as owner/
    // neighbour; after changeMesh that is a (owner, nei -1, region -1) boundary
    // face = the getFaceOrder FATAL. Remap any removed chain cell to its column's
    // remainder (a valid cell) so the residual face is merely poor (open/skew) and
    // the quality-retry orchestrator demotes it, instead of crashing the run.
    // No-op when the owner is an unsplit cell (cellCol miss) => external aero (no
    // dangling lateral faces, badFaces 0) is byte-identical.
    auto latRemap = [&](const label celli) -> label
    {
        if (celli >= 0 && cellCol.found(celli))
        {
            const label oc = cellCol[celli];
            if (!columns[oc].skip) { return stackOf[oc][columns[oc].nCol]; }
        }
        return celli;
    };
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
        const label rawOwn = mesh.faceOwner()[fi];
        const label rawNei =
            mesh.isInternalFace(fi) ? mesh.faceNeighbour()[fi] : -1;
        const label ownO = latRemap(rawOwn);
        const label neiO = latRemap(rawNei);
        // task 19: if a side was a REMOVED cell we just remapped (a transition
        // victim — the residual gap/skip-boundary faces), ATTRIBUTE the face to
        // that column so the quality-retry orchestrator can DEMOTE it. A lateral
        // face is otherwise column -1 = un-demotable, which left the manifold
        // "stuck" (split-owned failures -> 0, but these lateral ones persisted).
        // Untouched lateral faces (unsplit owners) stay column -1, unchanged.
        // recordFace only writes the debug map — the mesh is identical either way.
        label latCol = -1;
        if (rawOwn >= 0 && cellCol.found(rawOwn)) { latCol = cellCol[rawOwn]; }
        else if (rawNei >= 0 && cellCol.found(rawNei)) { latCol = cellCol[rawNei]; }
        const face nfF{labelList(nf)};
        saneCheck(nfF, fi, -1, "lateral");
        const label latPatch =
            mesh.isInternalFace(fi) ? -1 : patches.whichPatch(fi);
        badFace("lateral", fi, latCol, ownO, neiO, latPatch, -1, false);
        meshMod.modifyFace
        (
            nfF, fi, ownO, neiO, false, latPatch, -1, false
        );
        faceDone.set(fi);
        recordFace(fi, latCol >= 0 ? "lateralRemap" : "lateral", latCol);
    }

    // post-pass: any ORIGINAL face of a removed chain cell that no pass handled
    // (not in faceDone) dangles — it keeps its removed-cell owner and, if its
    // other side is also removed, ends as (owner, neighbour -1, region -1): the
    // changeMesh boundary-face crash. This catches the untouched-face case the
    // per-site badFace check cannot see.
    label nUnhandled = 0;
    forAllConstIters(cellCol, itc)
    {
        for (const label fi : mesh.cells()[itc.key()])
        {
            if (!faceDone.test(fi))
            {
                ++nUnhandled;
                if (nUnhandled <= 16)
                {
                    Info<< "SPLITLAYERS|UNHANDLED|fi|" << fi
                        << "|cell|" << itc.key() << "|ci|" << itc.val()
                        << "|internalFi|" << mesh.isInternalFace(fi)
                        << "|whichPatch|" << patches.whichPatch(fi) << nl;
                }
            }
        }
    }

    // ---- apply ------------------------------------------------------------
    Info<< "SPLITLAYERS|OOBtotal|" << nOOB << "|badFaces|" << nBad
        << "|deferredTop|" << nDefer << "|unhandled|" << nUnhandled << nl;
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
    mesh.updateMesh(map());

    // I-INT pinpoint: open cells (face-area-vectors don't sum to ~0) are the
    // variable-n termination defect. Attribute each to the constructs that built
    // its faces (final->origin faceMap + recordFace tags) + the owning column and
    // its nCol — so we see the deep/shallow boundary that fails to close.
    {
        const labelList& fMap = map().faceMap();
        const vectorField& fa = mesh.faceAreas();
        const labelList& fo = mesh.faceOwner();
        label nOpenCell = 0;
        // task 19: dump a column's seg0 ring-point coordinates (from newPtPos,
        // which changeMesh does not move; nOptSweeps == 0) per rail, with each
        // rail's wall base point. Guards skip / missing railKey / missing point
        // so the diagnostic can never throw.
        auto dumpSeg0 = [&](const label X, const char* side)
        {
            if (X < 0 || X >= columns.size() || columns[X].skip) { return; }
            const Column& cx = columns[X];
            const bool sharer = sharerRings.found(X);
            if (!sharer && !railsOf.found(railKey(cx, 0))) { return; }
            const labelList& asg = assignFor(X, cx);
            forAll(cx.rails, rc)
            {
                if (!sharer && !railsOf.found(railKey(cx, rc))) { continue; }
                const labelList& rg = ringsFor(X, cx, rc);
                const point& wp = pts[cx.rails[rc][0]];
                Info<< "SPLITLAYERS|OPPWALL|" << side << "|" << X
                    << "|rail|" << rc
                    << "|wall|(" << wp.x() << ' ' << wp.y() << ' '
                    << wp.z() << ")|seg0";
                forAll(asg, k)
                {
                    if (asg[k] != 0 || k >= rg.size()) { continue; }
                    const auto pit = newPtPos.cfind(rg[k]);
                    if (pit.good())
                    {
                        const point& rp = pit.val();
                        Info<< "|" << k << ":(" << rp.x() << ' '
                            << rp.y() << ' ' << rp.z() << ')';
                    }
                }
                Info<< nl;
            }
        };
        forAll(mesh.cells(), celli)
        {
            const cell& c = mesh.cells()[celli];
            vector s = Zero;
            scalar amax = VSMALL;
            for (const label fi : c)
            {
                s += (fo[fi] == celli ? 1.0 : -1.0) * fa[fi];
                amax = max(amax, mag(fa[fi]));
            }
            if (mag(s) > 1e-6 * amax)
            {
                ++nOpenCell;
                if (nOpenCell <= 16)
                {
                    word tags; label col = -1;
                    for (const label fi : c)
                    {
                        const label of = (fi < fMap.size()) ? fMap[fi] : -1;
                        if (of >= 0 && faceTagOf.found(of))
                        {
                            tags = tags + faceTagOf[of] + "|";
                            if (col < 0 && faceColOf.found(of))
                            {
                                col = faceColOf[of];
                            }
                        }
                    }
                    Info<< "SPLITLAYERS|OPENCELL|openness|" << mag(s)/amax
                        << "|col|" << col
                        << "|nCol|" << (col >= 0 ? columns[col].nCol : -1)
                        << "|tags|" << tags << nl;

                    // task 19: for an opposed open cell, dump the GEOMETRY of the
                    // two facing-wall columns' seg0 ring points so the collision is
                    // diagnosed from coordinates, not guessed. Find a strip face in
                    // this cell that names an opposing column (faceOtherColOf).
                    if (debugDump)
                    {
                        label colA = -1, colB = -1;
                        for (const label fi : c)
                        {
                            const label of = (fi < fMap.size()) ? fMap[fi] : -1;
                            if (of >= 0 && faceOtherColOf.found(of)
                             && faceColOf.found(of))
                            {
                                colA = faceColOf[of];
                                colB = faceOtherColOf[of];
                                break;
                            }
                        }
                        if (colA >= 0 && colB >= 0)
                        {
                            Info<< "SPLITLAYERS|OPPWALL|cell|" << celli
                                << "|openness|" << mag(s)/amax
                                << "|colA|" << colA
                                << "|nA|" << columns[colA].nCol
                                << "|colB|" << colB
                                << "|nB|" << columns[colB].nCol << nl;
                            dumpSeg0(colA, "A");
                            dumpSeg0(colB, "B");
                        }
                    }
                }
            }
        }
        Info<< "SPLITLAYERS|openCells|" << nOpenCell << nl;
    }

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

        // I3.1 quality-objective gate (accept-only-if-better). The smoothing
        // CANDIDATE is unchanged; each move is taken only if NEITHER the local max
        // nonOrtho NOR the local max skew rises. faceQ recomputes a face's nonOrtho
        // + skewness from a trial points field using checkMesh's own definitions
        // (cell volume-centroids via cell::centre, face geom via areaNormal/centre).
        // Per-point local-max non-increasing => global max non-increasing for BOTH
        // metrics (never-worse) — the property the old blind smoother lacked
        // (it regressed bracket 43.9->78.9, why sweeps shipped 0).
        const faceList& mFaces = mesh.faces();
        const cellList& mCells = mesh.cells();
        const labelList& fOwn = mesh.faceOwner();
        const labelList& fNei = mesh.faceNeighbour();
        const label nInt = mesh.nInternalFaces();

        auto faceQ = [&](const label fi, const pointField& P) -> Pair<scalar>
        {
            const point Co = mCells[fOwn[fi]].centre(P, mFaces);
            const point Cn = mCells[fNei[fi]].centre(P, mFaces);
            const vector Sf = mFaces[fi].areaNormal(P);
            const point Cf = mFaces[fi].centre(P);
            const vector d = Cn - Co;
            const scalar magd = mag(d), magS = mag(Sf);
            scalar no = 0, sk = 0;
            if (magd > SMALL && magS > SMALL)
            {
                no = radToDeg(Foam::acos(min(scalar(1), mag((d & Sf)/(magd*magS)))));
                const scalar den = (d & Sf);
                if (mag(den) > SMALL)
                {
                    const point Ci = Co + (((Cf - Co) & Sf)/den) * d;
                    sk = mag(Cf - Ci) / magd;
                }
            }
            return Pair<scalar>(no, sk);
        };
        // The affected set when np moves = ALL faces of the cells incident to np:
        // moving np shifts those cells' volume-centroids, which changes the nonOrtho
        // of EVERY one of their faces (owner−neighbour vector), not only the faces
        // that touch np. Gating on pointFaces alone missed the cells' opposite faces
        // and let the global max rise — so evaluate over pointCells' faces.
        auto localJ = [&](const label np, const pointField& P) -> Pair<scalar>
        {
            scalar wno = 0, wsk = 0;
            for (const label ci : mesh.pointCells()[np])
            {
                for (const label fi : mCells[ci])
                {
                    if (fi >= nInt) { continue; }
                    const Pair<scalar> q = faceQ(fi, P);
                    wno = max(wno, q.first()); wsk = max(wsk, q.second());
                }
            }
            return Pair<scalar>(wno, wsk);
        };
        auto maxLayerNonOrtho = [&]() -> scalar
        {
            scalar m = 0;
            forAll(inserted, i)
            {
                for (const label ci : mesh.pointCells()[inserted[i]])
                {
                    for (const label fi : mCells[ci])
                    {
                        if (fi < nInt) { m = max(m, faceQ(fi, newPts).first()); }
                    }
                }
            }
            return m;
        };

        const scalar noBefore = maxLayerNonOrtho();
        label nMoved = 0;
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
                const point bp = newPts[baseNewOf[np]];
                const point oldPos = newPts[np];
                const vector vOwn = oldPos - bp;
                const scalar L = mag(vOwn);
                if (L < SMALL) { continue; }

                // Local descent on the quality objective: try a small set of offset
                // DIRECTIONS — the Laplacian smooth + axis-tilt probes — each
                // renormalized to |vOwn| so the layer height (y+) is preserved, plus
                // "no move". Pick the candidate with the lowest nonOrtho that worsens
                // NEITHER nonOrtho nor skew vs the current position (so every move is
                // an improvement or a hold = never-worse). This actively reduces the
                // worst face instead of gating one fixed Laplacian candidate.
                const Pair<scalar> J0 = localJ(np, newPts);
                DynamicList<vector> dirs(8);
                dirs.append(vOwn);                                  // no move (never-worse floor)
                dirs.append(0.5 * vOwn + 0.5 * acc[np] / cnt[np]);  // Laplacian smooth
                const scalar dStep = 0.2 * L;
                for (direction a = 0; a < 3; ++a)
                {
                    if (geomD[a] < 0) { continue; }
                    vector e(Zero); e[a] = dStep;
                    dirs.append(vOwn + e);
                    dirs.append(vOwn - e);
                }
                point best = oldPos;
                Pair<scalar> bestJ = J0;
                for (const vector& vd : dirs)
                {
                    const scalar m = mag(vd);
                    if (m < SMALL) { continue; }
                    point cand = bp + vd * (L / m);                 // renormalize -> height kept
                    for (direction a = 0; a < 3; ++a)
                    {
                        if (geomD[a] < 0) { cand[a] = oldPos[a]; }  // freeze out-of-plane
                    }
                    newPts[np] = cand;
                    const Pair<scalar> Jc = localJ(np, newPts);
                    if (Jc.first()  <= J0.first()  + 1e-9           // never-worse: nonOrtho
                     && Jc.second() <= J0.second() + 1e-9           // never-worse: skew
                     && Jc.first()  <  bestJ.first() - 1e-12)       // strictly better nonOrtho
                    {
                        best = cand; bestJ = Jc;
                    }
                }
                newPts[np] = best;
                if (mag(best - oldPos) > SMALL) { ++nMoved; }
            }
        }
        const scalar noAfter = maxLayerNonOrtho();

        // WORSTFACE locator (behaviour-neutral): the highest-nonOrtho layer faces,
        // each with its construct tag (built->orig faceMap + faceTagOf, populated
        // under -debugDump) + centre — to confirm whether the standing max sits at
        // sharp EDGES (ridge/strip constructs) before any edge-aware construction
        // (v3). At sweeps>0 these are the faces the descent could NOT reduce.
        {
            const labelList& fMap = map().faceMap();
            DynamicList<scalar> noList;
            DynamicList<label> fidList;
            labelHashSet seen;
            forAll(inserted, i)
            {
                for (const label ci : mesh.pointCells()[inserted[i]])
                {
                    for (const label fi : mCells[ci])
                    {
                        if (fi < nInt && seen.insert(fi))
                        {
                            noList.append(faceQ(fi, newPts).first());
                            fidList.append(fi);
                        }
                    }
                }
            }
            labelList order;
            sortedOrder(noList, order);                       // ascending
            const label n = order.size();
            const label nShow = min(label(20), n);
            for (label k = 0; k < nShow; ++k)
            {
                const label idx = order[n - 1 - k];           // worst first
                const label fi = fidList[idx];
                const point Cf = mesh.faces()[fi].centre(newPts);
                word tag("none"); label col = -1;
                const label of = (fi < fMap.size() ? fMap[fi] : -1);
                if (of >= 0 && faceTagOf.found(of))
                {
                    tag = faceTagOf[of];
                    if (faceColOf.found(of)) { col = faceColOf[of]; }
                }
                Info<< "SPLITLAYERS|WORSTFACE|rank|" << k
                    << "|nonOrtho|" << noList[idx]
                    << "|centre|(" << Cf.x() << ' ' << Cf.y() << ' ' << Cf.z() << ")"
                    << "|tag|" << tag << "|col|" << col << nl;
            }
        }

        mesh.movePoints(newPts);
        Info<< "SPLITLAYERS|optimize|sweeps|" << nOptSweeps
            << "|points|" << inserted.size() << "|moved|" << nMoved
            << "|nonOrthoBefore|" << noBefore << "|nonOrthoAfter|" << noAfter << nl;
    }

    mesh.setInstance(runTime.constant());
    mesh.write();

    Info<< "SPLITLAYERS|result|split-ok" << nl << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
