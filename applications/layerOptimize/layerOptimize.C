/*---------------------------------------------------------------------------*\
  layerOptimize — unified never-worse layer post-smoother (H2b).

  Design: mesh-tools-gpl/docs/2026-06-21-layeroptimize-unified-postsmoother.md

  Gives the EXTRUSION arm (and any layer polyMesh) the same proven optimizer the
  split arm carries inline. Reads a polyMesh + a layerCells cellSet + the wall
  patch; recovers the prism columns from topology (walk each wall face's stack
  outward via opposite faces, tracking each layer point's wall base via the hex
  vertical edges); then runs the ported localDescent — active-set, faceQ
  never-worse descent (nonOrtho + skew monotonically non-increasing), y+-pinned
  (each move renormalized to keep the wall-offset magnitude), early-stopped.

  v1: localDescent only (cornerScale is a follow-on). -detectOnly = the spike
  (column detection + stats, no move). faceQ/localJ/the descent are copied from
  splitLayers.C; extract to a shared header once green.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "unitConversion.H"
#include "Pair.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption("detectOnly", "column detection + stats only (no move)");
    argList::addOption("wallPatch", "name", "BL wall patch (default: largest type-wall patch)");
    argList::addOption("cellSet", "name", "layer-cell set (default: layerCells)");
    argList::addOption("optimizeSweeps", "int", "max descent sweeps (default 10)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word cellSetName(args.getOrDefault<word>("cellSet", "layerCells"));
    const label nOptSweeps = args.getOrDefault<label>("optimizeSweeps", 10);

    // ── resolve the wall patch (explicit, else largest type-wall patch) ──────
    label patchID = -1;
    if (args.found("wallPatch"))
    {
        patchID = mesh.boundaryMesh().findPatchID(args.get<word>("wallPatch"));
    }
    else
    {
        label bestSz = -1;
        forAll(mesh.boundaryMesh(), pi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[pi];
            if (pp.type() == "wall" && pp.size() > bestSz)
            {
                bestSz = pp.size(); patchID = pi;
            }
        }
    }
    if (patchID < 0)
    {
        FatalErrorInFunction << "no wall patch (pass -wallPatch). patches = "
            << mesh.boundaryMesh().names() << exit(FatalError);
    }
    const polyPatch& wallP = mesh.boundaryMesh()[patchID];
    const cellSet layerCells(mesh, cellSetName);

    const faceList& mFaces = mesh.faces();
    const cellList& mCells = mesh.cells();
    const labelList& fOwn = mesh.faceOwner();
    const labelList& fNei = mesh.faceNeighbour();
    const label nInt = mesh.nInternalFaces();

    auto oppositeFace = [&](const label c, const label fIn) -> label
    {
        labelHashSet inPts;
        for (const label p : mFaces[fIn]) { inPts.insert(p); }
        for (const label f : mCells[c])
        {
            if (f == fIn) { continue; }
            bool shares = false;
            for (const label p : mFaces[f]) { if (inPts.found(p)) { shares = true; break; } }
            if (!shares) { return f; }
        }
        return -1;
    };
    // vertical-edge partner of fOut-point p in the inner face = the inner-face
    // point adjacent to p across a side face (the hex's vertical edge).
    auto wallPartner = [&](const label c, const labelHashSet& inPts, const label p) -> label
    {
        for (const label f : mCells[c])
        {
            const face& fc = mFaces[f];
            forAll(fc, j)
            {
                if (fc[j] == p)
                {
                    const label a = fc[fc.fcIndex(j)];
                    const label b = fc[fc.rcIndex(j)];
                    if (inPts.found(a)) { return a; }
                    if (inPts.found(b)) { return b; }
                }
            }
        }
        return -1;
    };

    // ── walk the prism columns: build wallBaseOf + depthOf ───────────────────
    Map<label> wallBaseOf;   // point -> its wall base point (y+ origin)
    Map<label> depthOf;      // point -> layer index k (0 = wall, fixed)
    label nCols = 0, nNonHex = 0, nNoWall = 0;
    DynamicList<label> depths;

    forAll(wallP, i)
    {
        const label wf = wallP.start() + i;
        label c = fOwn[wf];
        if (!layerCells.found(c)) { ++nNoWall; continue; }
        labelHashSet inPts;
        for (const label p : mFaces[wf])
        {
            inPts.insert(p);
            if (!depthOf.found(p)) { wallBaseOf.set(p, p); depthOf.set(p, 0); }
        }
        label fIn = wf, k = 0;
        bool nonHex = false;
        while (true)
        {
            if (mCells[c].size() != 6) { nonHex = true; break; }
            const label fOut = oppositeFace(c, fIn);
            if (fOut < 0) { break; }
            ++k;
            for (const label p : mFaces[fOut])
            {
                if (!depthOf.found(p) || k < depthOf[p])   // shortest (most-vertical) column wins
                {
                    const label q = wallPartner(c, inPts, p);
                    if (q >= 0 && wallBaseOf.found(q))
                    {
                        wallBaseOf.set(p, wallBaseOf[q]);
                        depthOf.set(p, k);
                    }
                }
            }
            if (fOut >= nInt) { break; }
            const label cn = (fOwn[fOut] == c) ? fNei[fOut] : fOwn[fOut];
            if (!layerCells.found(cn)) { break; }
            inPts.clear();
            for (const label p : mFaces[fOut]) { inPts.insert(p); }
            c = cn; fIn = fOut;
        }
        if (nonHex) { ++nNonHex; continue; }
        ++nCols; depths.append(k);
    }

    DynamicList<label> movable;          // layer points (k>=1); wall points (k=0) stay fixed
    forAllConstIters(depthOf, it) { if (it.val() >= 1) { movable.append(it.key()); } }

    label dmin = labelMax, dmax = 0; scalar dsum = 0;
    forAll(depths, i) { dmin = min(dmin, depths[i]); dmax = max(dmax, depths[i]); dsum += depths[i]; }
    Info<< "LAYEROPT|detect|wallPatch|" << wallP.name() << "|wallFaces|" << wallP.size()
        << "|columns|" << nCols << "|nonHexSkipped|" << nNonHex << "|noLayerCellAtWall|" << nNoWall
        << "|movablePts|" << movable.size() << "|depthMin|" << (nCols ? dmin : 0)
        << "|depthMax|" << dmax << "|depthMean|" << (nCols ? dsum/nCols : 0.0) << nl;

    if (args.found("detectOnly") || nOptSweeps <= 0 || movable.empty())
    {
        Info<< nl << "End" << nl;
        return 0;
    }

    // ── ported quality objective (faceQ/localJ) — checkMesh definitions ───────
    pointField newPts(mesh.points());
    const labelListList& pointPoints = mesh.pointPoints();
    const Vector<label>& geomD = mesh.geometricD();

    auto faceQ = [&](const label fi, const pointField& P) -> Pair<scalar>
    {
        const point Co = mCells[fOwn[fi]].centre(P, mFaces);
        const point Cn = mCells[fNei[fi]].centre(P, mFaces);
        const vector Sf = mFaces[fi].areaNormal(P);
        const point Cf = mFaces[fi].centre(P);
        const vector d = Cn - Co; const scalar magd = mag(d), magS = mag(Sf);
        scalar no = 0, sk = 0;
        if (magd > SMALL && magS > SMALL)
        {
            // checkMesh definition, SIGNED (0..180): no mag() — folding 90..180 to
            // 90..0 would make the gate REWARD inversion (the v1 tangle bug). With
            // the true angle, any move past 90 is correctly seen as worse → rejected.
            no = radToDeg(Foam::acos(max(scalar(-1), min(scalar(1), (d & Sf)/(magd*magS)))));
            const scalar den = (d & Sf);
            if (mag(den) > SMALL)
            { const point Ci = Co + (((Cf - Co) & Sf)/den) * d; sk = mag(Cf - Ci)/magd; }
        }
        return Pair<scalar>(no, sk);
    };
    auto localJ = [&](const label np, const pointField& P) -> Pair<scalar>
    {
        scalar wno = 0, wsk = 0;
        for (const label ci : mesh.pointCells()[np])
            for (const label fi : mCells[ci])
                if (fi < nInt)
                { const Pair<scalar> q = faceQ(fi, P); wno = max(wno, q.first()); wsk = max(wsk, q.second()); }
        return Pair<scalar>(wno, wsk);
    };
    auto maxLayerNonOrtho = [&]() -> scalar
    {
        scalar m = 0;
        forAll(movable, i)
            for (const label ci : mesh.pointCells()[movable[i]])
                for (const label fi : mCells[ci])
                    if (fi < nInt) { m = max(m, faceQ(fi, newPts).first()); }
        return m;
    };

    const scalar OPT_ACTIVE_FRAC = 0.7, OPT_ACTIVE_FLOOR = 30.0, OPT_CONV_TOL = 1e-4;
    const scalar noBefore = maxLayerNonOrtho();

    label nMoved = 0; scalar prevMax = VGREAT;
    for (label sweep = 0; sweep < nOptSweeps; ++sweep)
    {
        // Laplacian accumulator over same-layer neighbours (y+-band preserved)
        vectorField acc(newPts.size(), Zero);
        scalarField cnt(newPts.size(), 0.0);
        forAll(movable, i)
        {
            const label np = movable[i];
            for (const label nb : pointPoints[np])
                if (depthOf.found(nb) && depthOf[nb] == depthOf[np])
                { acc[np] += newPts[nb] - newPts[wallBaseOf[nb]]; cnt[np] += 1.0; }
        }

        // active set: movable points incident to a >=T face (T tracks curMax)
        boolList activePt(newPts.size(), false);
        scalar curMax = 0;
        {
            DynamicList<label> lf; DynamicList<scalar> lno; labelHashSet seenA;
            forAll(movable, i)
                for (const label ci : mesh.pointCells()[movable[i]])
                    for (const label fi : mCells[ci])
                        if (fi < nInt && seenA.insert(fi))
                        { const scalar no = faceQ(fi, newPts).first(); lf.append(fi); lno.append(no); curMax = max(curMax, no); }
            const scalar T = max(OPT_ACTIVE_FLOOR, OPT_ACTIVE_FRAC*curMax);
            forAll(lf, j)
                if (lno[j] >= T)
                {
                    const label fi = lf[j];
                    for (const label vp : mCells[fOwn[fi]].labels(mFaces)) { activePt[vp] = true; }
                    for (const label vp : mCells[fNei[fi]].labels(mFaces)) { activePt[vp] = true; }
                }
        }
        if (sweep > 0 && (prevMax - curMax) < OPT_CONV_TOL*max(prevMax, scalar(SMALL)))
        {
            Info<< "LAYEROPT|converged|sweep|" << sweep << "|curMax|" << curMax << nl;
            break;
        }
        prevMax = curMax;

        forAll(movable, i)
        {
            const label np = movable[i];
            if (!activePt[np] || cnt[np] < 1) { continue; }
            const point bp = newPts[wallBaseOf[np]];
            const point oldPos = newPts[np];
            const vector vOwn = oldPos - bp; const scalar L = mag(vOwn);
            if (L < SMALL) { continue; }
            const Pair<scalar> J0 = localJ(np, newPts);
            DynamicList<vector> dirs(8);
            dirs.append(vOwn);
            dirs.append(0.5*vOwn + 0.5*acc[np]/cnt[np]);
            const scalar dStep = 0.2*L;
            for (direction a = 0; a < 3; ++a)
            {
                if (geomD[a] < 0) { continue; }
                vector e(Zero); e[a] = dStep;
                dirs.append(vOwn + e); dirs.append(vOwn - e);
            }
            point best = oldPos; Pair<scalar> bestJ = J0;
            for (const vector& vd : dirs)
            {
                const scalar m = mag(vd);
                if (m < SMALL) { continue; }
                point cand = bp + vd*(L/m);                    // renormalize -> y+ kept
                for (direction a = 0; a < 3; ++a) { if (geomD[a] < 0) { cand[a] = oldPos[a]; } }
                newPts[np] = cand;
                const Pair<scalar> Jc = localJ(np, newPts);
                if (Jc.first() <= J0.first() + 1e-9 && Jc.second() <= J0.second() + 1e-9
                 && Jc.first() < bestJ.first() - 1e-12)
                { best = cand; bestJ = Jc; }
            }
            newPts[np] = best;
            if (mag(best - oldPos) > SMALL) { ++nMoved; }
        }
    }

    const scalar noAfter = maxLayerNonOrtho();
    mesh.movePoints(newPts);
    mesh.write();
    Info<< "LAYEROPT|optimize|sweeps|" << nOptSweeps << "|movable|" << movable.size()
        << "|moved|" << nMoved << "|nonOrthoBefore|" << noBefore
        << "|nonOrthoAfter|" << noAfter << nl;

    Info<< nl << "End" << nl;
    return 0;
}


// ************************************************************************* //
