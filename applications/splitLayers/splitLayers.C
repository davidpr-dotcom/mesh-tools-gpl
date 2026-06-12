/*---------------------------------------------------------------------------*\
    splitLayers — split-based boundary-layer generation (WP 5.4-spike)

    Carves graded boundary layers out of existing wall-adjacent hex cells by
    inserting points along the wall-normal cell columns and rebuilding the
    topology. No mesh movement of pre-existing points, no extrusion, no
    coverage failure mode: every wall face of the target patch that heads a
    hex column receives exactly the specified layers (columns that cannot be
    split are skipped and REPORTED — the skip rate is a primary spike metric).

    Pipeline:
      1. spec        layer-spec JSON (etc/contracts/layer-spec, semver 0.x)
      2. columns     wall face -> opposing face of the owner hex; the four
                     point-to-point cell edges form the column
      3. placement   ring points at absolute heights y1*er^k along either
                     the column edge (-placement edge) or the smoothed
                     inward point-normal field (-placement normal, default;
                     falls back to the edge direction beyond 60 degrees)
      4. topology    polyTopoChange: n layer cells + remainder per column;
                     side faces of split cells become strip stacks; lateral
                     faces of non-split neighbours keep one face with the
                     ring points inserted into their polygon (no hanging
                     points)
      5. optimize    -optimizeSweeps N (default 10, 0 = off): tangential
                     Laplacian smoothing of ring points with the radial
                     height constraint re-imposed each sweep (wall points
                     fixed by construction)

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
#include "processorPolyPatch.H"
#include "EdgeMap.H"
#include "edgeList.H"
#include "IFstream.H"
#include "unitConversion.H"

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

// Area-weighted inward point normals on the patch, Laplacian-smoothed over
// patch points. Feature points (adjacent patch-face normals differing by
// more than featureAngle) take part like any other point — smoothing is
// what produces the bisector-like directions there.
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

    // feature-point count (diagnostic): max angle between adjacent faces
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

    // Laplacian smoothing over patch point-point connectivity
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

// * * * * * * * * * * * * * column structure  * * * * * * * * * * * * * * //

struct Column
{
    label patchFacei;       // local patch face index
    label meshFacei;        // global wall face
    label celli;            // owner hex
    label oppFacei;         // opposing face
    FixedList<label, 4> base;   // wall-face points (mesh labels)
    FixedList<label, 4> top;    // matching opposing points (mesh labels)
};

int main(int argc, char *argv[])
{
    argList::addNote("Split-based boundary layers (WP 5.4-spike, full algorithm)");
    argList::addOption("spec", "file", "layer-spec JSON (contract 0.x)");
    argList::addOption("placement", "edge|normal",
                       "ring-point placement (default: normal)");
    argList::addOption("optimizeSweeps", "N",
                       "tangential smoothing sweeps (default 10, 0=off)");
    argList::addOption("featureAngle", "deg",
                       "feature-point detection angle (default 45)");
    argList::addBoolOption("dryRun", "identify and report columns only");
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Create mesh for time = " << runTime.timeName() << nl << endl;
    // polyMesh (not fvMesh): a pure mesh tool must not create solver-side
    // artifacts (fvMesh::movePoints writes 0/meshPhi, which then breaks
    // any subsequent tool run on the case — found 2026-06-12).
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

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelListList& cellEdges = mesh.cellEdges();
    const edgeList& meshEdges = mesh.edges();

    // mark cells owned by >1 target patch face (corner cells) -> skip
    labelList patchFaceCount(mesh.nCells(), 0);

    // ---- gather per-patch columns --------------------------------------
    List<DynamicList<Column>> allColumns(specs.size());
    List<label> skippedNonHex(specs.size(), 0);
    List<label> skippedCorner(specs.size(), 0);

    forAll(specs, si)
    {
        const label patchi = patches.findPatchID(specs[si].patchName);
        if (patchi < 0)
        {
            FatalErrorInFunction
                << "patch '" << specs[si].patchName << "' not found"
                << exit(FatalError);
        }
        const polyPatch& pp = patches[patchi];
        for (const label c : pp.faceCells()) { patchFaceCount[c]++; }
    }

    forAll(specs, si)
    {
        const polyPatch& pp = patches[patches.findPatchID(specs[si].patchName)];

        forAll(pp, pfi)
        {
            const label fG = pp.start() + pfi;
            const label own = mesh.faceOwner()[fG];

            if (patchFaceCount[own] > 1)
            {
                ++skippedCorner[si];
                continue;
            }
            const cell& c = cells[own];
            if (c.size() != 6 || faces[fG].size() != 4)
            {
                ++skippedNonHex[si];
                continue;
            }
            const label oppFi = c.opposingFaceLabel(fG, faces);
            if (oppFi < 0 || faces[oppFi].size() != 4)
            {
                ++skippedNonHex[si];
                continue;
            }

            Column col;
            col.patchFacei = pfi;
            col.meshFacei = fG;
            col.celli = own;
            col.oppFacei = oppFi;

            const face& wf = faces[fG];
            const face& of = faces[oppFi];
            bool ok = true;
            forAll(wf, k)
            {
                col.base[k] = wf[k];
                // column edge: cell edge from wf[k] to a point of opp face
                label qFound = -1;
                for (const label ei : cellEdges[own])
                {
                    const edge& e = meshEdges[ei];
                    const label other = e.otherVertex(wf[k]);
                    if (other >= 0 && of.found(other))
                    {
                        qFound = other;
                        break;
                    }
                }
                if (qFound < 0) { ok = false; break; }
                col.top[k] = qFound;
            }
            if (!ok)
            {
                ++skippedNonHex[si];
                continue;
            }
            allColumns[si].append(col);
        }

        Info<< "SPLITLAYERS|patch|" << specs[si].patchName
            << "|nFaces|" << pp.size()
            << "|columns|" << allColumns[si].size()
            << "|skippedNonHex|" << skippedNonHex[si]
            << "|skippedCorner|" << skippedCorner[si]
            << nl;
    }

    if (dryRun)
    {
        Info<< "SPLITLAYERS|result|dry-run-ok" << nl << "End\n" << endl;
        return 0;
    }

    // ---- ring point creation (shared per column edge) -------------------
    polyTopoChange meshMod(mesh);

    // per spec: heights t_k
    List<scalarList> cumHeights(specs.size());
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
    }

    // smoothed normals per patch (point field, patch-local indexing)
    List<vectorField> patchNormals(specs.size());
    forAll(specs, si)
    {
        const polyPatch& pp = patches[patches.findPatchID(specs[si].patchName)];
        label nFeat = 0;
        patchNormals[si] = smoothedInwardNormals(pp, 5, featAngle, nFeat);
        Info<< "SPLITLAYERS|normals|" << specs[si].patchName
            << "|featurePoints|" << nFeat << nl;
    }

    // EdgeMap: column edge (p_base, p_top) -> list of new ring point labels
    EdgeMap<labelList> ringPoints;
    EdgeMap<label> edgeSpec;               // column edge -> spec index

    const pointField& pts = mesh.points();

    forAll(specs, si)
    {
        const polyPatch& pp = patches[patches.findPatchID(specs[si].patchName)];
        const scalar cosGuard = Foam::cos(degToRad(60.0));

        for (const Column& col : allColumns[si])
        {
            forAll(col.base, k)
            {
                const edge colEdge(col.base[k], col.top[k]);
                if (ringPoints.found(colEdge)) { continue; }

                const point& p = pts[col.base[k]];
                const point& q = pts[col.top[k]];
                const vector edgeDir = (q - p) / max(mag(q - p), SMALL);
                const scalar edgeLen = mag(q - p);

                vector dir = edgeDir;
                if (placement == "normal")
                {
                    const label loc = pp.whichPoint(col.base[k]);
                    if (loc >= 0)
                    {
                        const vector& ns = patchNormals[si][loc];
                        if ((ns & edgeDir) > cosGuard) { dir = ns; }
                    }
                }

                // compress if the column is too short for the stack
                const scalar total = cumHeights[si].last();
                const scalar scaleH =
                    (total < 0.9 * edgeLen) ? 1.0 : 0.9 * edgeLen / total;

                labelList ring(specs[si].nLayers);
                forAll(ring, j)
                {
                    const point rp = p + dir * (scaleH * cumHeights[si][j]);
                    ring[j] = meshMod.addPoint(rp, col.base[k], -1, true);
                }
                ringPoints.insert(colEdge, ring);
                edgeSpec.insert(colEdge, si);
            }
        }
    }

    // ---- cells: replace each split cell by n layers + remainder ---------
    // newCells[celli] = list of n+1 new cell labels (wall-side first)
    Map<labelList> newCells;

    forAll(specs, si)
    {
        for (const Column& col : allColumns[si])
        {
            labelList stack(specs[si].nLayers + 1);
            forAll(stack, j)
            {
                stack[j] = meshMod.addCell(-1, -1, -1, col.celli, -1);
            }
            meshMod.removeCell(col.celli, -1);
            newCells.insert(col.celli, stack);
        }
    }

    // helper: ring point for (column edge, layer k)
    auto ringAt = [&](const label basePt, const label topPt, const label j)
    {
        return ringPoints[edge(basePt, topPt)][j];
    };

    // ---- faces -----------------------------------------------------------
    // For each split cell: wall face -> layer cell 0 (modify);
    // ring faces between consecutive stack cells (add);
    // opposing face -> remainder (handled in the generic pass below);
    // side faces -> strips; lateral faces of unsplit neighbours -> point
    // insertion. The generic pass walks every face of every split cell once.

    bitSet faceDone(mesh.nFaces(), false);

    forAll(specs, si)
    {
        const label patchi = patches.findPatchID(specs[si].patchName);
        const label n = specs[si].nLayers;

        for (const Column& col : allColumns[si])
        {
            const labelList& stack = newCells[col.celli];
            const face& wf = faces[col.meshFacei];

            // 1. wall face -> owner = first layer cell
            if (!faceDone.test(col.meshFacei))
            {
                meshMod.modifyFace
                (
                    wf, col.meshFacei, stack[0], -1, false, patchi, -1, false
                );
                faceDone.set(col.meshFacei);
            }

            // 2. internal ring faces between stack cells
            for (label j = 0; j < n; ++j)
            {
                face rf(4);
                forAll(wf, k)
                {
                    // reversed order so the normal points away from the wall
                    rf[k] = ringAt(col.base[3 - k], col.top[3 - k], j);
                }
                meshMod.addFace
                (
                    rf, stack[j], stack[j + 1],
                    -1, -1, col.meshFacei, false, -1, -1, false
                );
            }

            // 3. remaining original faces of this cell
            const cell& c = cells[col.celli];
            for (const label fi : c)
            {
                if (fi == col.meshFacei || faceDone.test(fi)) { continue; }

                const face& f = faces[fi];
                const label ownO = mesh.faceOwner()[fi];
                const bool isOwner = (ownO == col.celli);
                const label other =
                    mesh.isInternalFace(fi)
                  ? (isOwner ? mesh.faceNeighbour()[fi] : ownO)
                  : -1;

                if (fi == col.oppFacei)
                {
                    // opposing face -> remainder cell
                    const label nei =
                        (other >= 0 && newCells.found(other))
                      ? newCells[other][specs[si].nLayers]  // unlikely; guard
                      : other;
                    if (isOwner)
                    {
                        meshMod.modifyFace(f, fi, stack[n], nei, false,
                            mesh.isInternalFace(fi) ? -1
                            : patches.whichPatch(fi), -1, false);
                    }
                    else
                    {
                        meshMod.modifyFace(f, fi, nei, stack[n], false,
                            -1, -1, false);
                    }
                    faceDone.set(fi);
                    continue;
                }

                // side face: contains exactly 2 base + 2 top points
                label ia = -1, ib = -1;
                forAll(col.base, k)
                {
                    if (f.found(col.base[k]))
                    {
                        (ia < 0 ? ia : ib) = k;
                    }
                }

                if (ia >= 0 && ib >= 0)
                {
                    // Strip-split this side face by POINT SUBSTITUTION on
                    // the original face: winding (and therefore the
                    // owner/neighbour orientation convention) is inherited
                    // from f. Levels: 0 = wall points, k = ring k-1,
                    // n+1 = top points; strip j spans level j..j+1.
                    const bool otherSplit =
                        (other >= 0 && newCells.found(other));

                    auto levelPoint =
                        [&](const label meshPt, const label lev) -> label
                    {
                        forAll(col.base, k)
                        {
                            if (meshPt == col.base[k])
                            {
                                return (lev == 0)
                                    ? col.base[k]
                                    : ringAt(col.base[k], col.top[k], lev-1);
                            }
                            if (meshPt == col.top[k])
                            {
                                return (lev == n + 1)
                                    ? col.top[k]
                                    : ringAt(col.base[k], col.top[k], lev-1);
                            }
                        }
                        return meshPt;   // not on a column edge (impossible
                                         // for a hex side face, but safe)
                    };

                    for (label j = 0; j <= n; ++j)
                    {
                        face sf(f.size());
                        forAll(f, k)
                        {
                            // wall-side points of f -> level j,
                            // top-side points of f -> level j+1
                            bool isBase = false;
                            forAll(col.base, m)
                            {
                                if (f[k] == col.base[m]) { isBase = true; }
                            }
                            sf[k] = levelPoint(f[k], isBase ? j : j + 1);
                        }

                        // ownership: replace celli by stack[j]; if the
                        // other side is split too, replace it by its j-th
                        // stack cell; otherwise keep the original label.
                        label sOwn = ownO;
                        label sNei = mesh.isInternalFace(fi)
                                   ? mesh.faceNeighbour()[fi] : -1;
                        if (sOwn == col.celli) { sOwn = stack[j]; }
                        else if (otherSplit && sOwn == other)
                        {
                            sOwn = newCells[other][j];
                        }
                        if (sNei == col.celli) { sNei = stack[j]; }
                        else if (otherSplit && sNei == other)
                        {
                            sNei = newCells[other][j];
                        }
                        const label sPatch = mesh.isInternalFace(fi)
                                           ? -1 : patches.whichPatch(fi);

                        if (j == 0)
                        {
                            meshMod.modifyFace(sf, fi, sOwn, sNei, false,
                                sPatch, -1, false);
                        }
                        else
                        {
                            meshMod.addFace(sf, sOwn, sNei,
                                -1, -1, fi, false, sPatch, -1, false);
                        }
                    }
                    faceDone.set(fi);
                }
            }
        }
    }

    // ---- lateral faces of UNSPLIT cells touching column edges -----------
    // Insert ring points into their polygons (no hanging nodes). A face may
    // touch SEVERAL split column edges (e.g. corner-cell faces between two
    // split neighbours) — collect the face set first, then walk each face
    // once and insert the rings of every column edge it contains.
    labelHashSet lateralFaces;
    forAllConstIters(ringPoints, iter)
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
            const auto it = ringPoints.cfind(fe);
            if (it.good())
            {
                const labelList& ring = it.val();
                if (f[k] == it.key()[0])    // stored as (base, top)
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
        meshMod.modifyFace
        (
            face(labelList(nf)), fi, ownO, neiO, false,
            mesh.isInternalFace(fi) ? -1 : patches.whichPatch(fi),
            -1, false
        );
        faceDone.set(fi);
    }

    // ---- apply -----------------------------------------------------------
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
    mesh.updateMesh(map());

    Info<< "SPLITLAYERS|topology|cells|" << mesh.nCells()
        << "|points|" << mesh.nPoints() << nl;

    // ---- node-position optimization (tangential Laplacian, radial fix) --
    if (nOptSweeps > 0)
    {
        pointField newPts(mesh.points());

        // Recover the ring points in FINAL numbering: every ring point was
        // added with its wall point as master, so pointMap (new->old) sends
        // it to the old wall-point label. Per column edge: collect all new
        // points mapping to the old base label, drop the base itself, sort
        // by distance from the (fixed) base point -> layer order.
        const labelList& pointMap = map().pointMap();
        const labelList& revPointMap = map().reversePointMap();

        Map<DynamicList<label>> children;   // old base label -> new points
        forAll(pointMap, np)
        {
            const label oldp = pointMap[np];
            if (oldp >= 0 && revPointMap[oldp] != np)
            {
                children(oldp).append(np);
            }
        }

        DynamicList<label> inserted;
        DynamicList<label> insBase, insK, insSpec;
        forAllConstIters(ringPoints, iter)
        {
            const label b = iter.key()[0];   // edge stores base point first
            if (!children.found(b)) { continue; }
            const label newBase = revPointMap[b];
            DynamicList<label>& kids = children[b];
            // sort by distance from base
            labelList sorted(kids);
            std::sort(sorted.begin(), sorted.end(),
                [&](label x, label y)
                {
                    return magSqr(newPts[x] - newPts[newBase])
                         < magSqr(newPts[y] - newPts[newBase]);
                });
            const label si = edgeSpec[iter.key()];
            forAll(sorted, j)
            {
                if (j >= specs[si].nLayers) { break; }
                inserted.append(sorted[j]);
                insBase.append(b);
                insK.append(j);
                insSpec.append(si);
            }
            children.erase(b);   // each base belongs to one column edge
        }

        // neighbour structure: ring points sharing base-point adjacency is
        // expensive to rebuild exactly; approximate neighbours = inserted
        // points within the same layer connected by a mesh edge
        const labelListList& pointPoints = mesh.pointPoints();
        Map<label> layerOf, specOf, baseOf;
        forAll(inserted, i)
        {
            layerOf.insert(inserted[i], insK[i]);
            specOf.insert(inserted[i], insSpec[i]);
            baseOf.insert(inserted[i], insBase[i]);
        }
        // Smooth the OFFSET VECTORS (point - own base), not positions:
        // neighbours exchange direction information only. Uniform regions
        // are a fixed point (no drift, no stack rotation); at corners and
        // feature edges directions blend — the intended behaviour.
        labelList baseNewOf(newPts.size(), -1);
        forAll(inserted, i)
        {
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
                if (mv > SMALL)
                {
                    // preserve this point's own (possibly compressed) height
                    v *= mag(vOwn) / mv;
                }
                else
                {
                    v = vOwn;
                }
                point target = bp + v;
                // never move along empty (non-geometric) directions (2D)
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
