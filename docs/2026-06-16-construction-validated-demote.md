# splitLayers: construction-validated demote (I1.1 convex-ridge fix + I1.2 rung-1, 2026-06-16)

> **SUPERSEDED 2026-06-16.** The in-tool detector this doc proposed does not work:
> an edge-aspect "sliver" test flags every legitimately thin boundary-layer cell
> (BL cells ARE high-aspect by design), and the real discriminator (checkMesh
> skewness/orientation) needs the built mesh, not construction-time geometry. The
> shipped fix moved validate→attribute→demote→retry to the **meshcore Python
> orchestrator** (`meshcore/executors/splitlayers_runner.py`) + a `splitLayers
> -skipColumns` option — driven by the real checkMesh result. It WORKS (cow: 79
> faces → 23 columns demoted → clean @ ~99.2%; bracket/ahmed recovered, clean).
> Kept for the record / the demote-and-rebuild reasoning. The §"demote, rebuild
> once" idea survives as the orchestrator's retry loop.

**Track:** Phase 5 / WP 5.2. Replaces the geometric `skippedConvexRidge` proxy
(I0.1 interim) with a **validate-the-real-geometry-then-demote** mechanism — the
first rung of the I1.2 fallback ladder, fused with the I4.1 self-check.
**Build/iterate on WSL** (`zefra-foam`, OF v2512); the rebuild wrap is the delicate
part — compile-iterate it.

## Why (the S0 + convex-ridge data, 2026-06-16)

Disabling the skip (`-convexRidgeAngle 179`) + checkMesh:
- **bracket** 1562→**1607 cols, Mesh OK**; **ahmed** 1707→**1891 cols, Mesh OK**
  (all 184 ridge-skips were FALSE POSITIVES); **cow** 21214→21425 but checkMesh
  **FAILED** (30 wrong-oriented + 40 skew — the I0.1 `strips:opposed` defect).

The geometric skip is the **wrong kind of check**: it pre-emptively guesses from
wall-edge angle, so it over-fires wherever a sharp convex edge exists (bracket/
ahmed, where the strips are actually clean) and a single global angle cannot
separate those from the cow ears (where the strips are genuinely degenerate). The
robust, polyvalent fix is to **stop guessing and validate the built result**: keep
a column only if the geometry it actually produces is sane.

## Design — validate the real faces, demote, rebuild once

The tool already has the detector: the construction-time `saneCheck` lambda (spread
test + the I4.1 sliver test) flags degenerate faces **with their producing column**
(`SPLITLAYERS|INSANE|...|col|<ci>`). Today it only *reports*. Make it *drive a
demote*:

1. **Collect, don't just report.** `saneCheck` inserts the offending column into a
   `labelHashSet failedCols` (it already has `ci2`).
2. **Single rebuild.** Wrap construction (PASS 3 ring-creation → PASS 4 topology →
   the lateral/absorption passes, everything from `polyTopoChange meshMod(mesh)`
   down to just before `mesh.setInstance/write`) in a **2-iteration loop**. After
   iteration 0, if `failedCols` is non-empty: set `col.skip=true` for each (counted
   as `nDemoted`), **reset** the per-construction state, and run iteration 1 with
   them demoted. Iteration 1 builds only sane columns.
3. **Backstop.** The I4.1 self-check (refuse-on-`nInsane>0`) stays: if anything is
   still INSANE after the single rebuild, refuse (`split-FAILED`) rather than write
   a bad mesh. (Two clean iterations is the expected path; the backstop covers the
   pathological residual.)
4. **Remove the geometric convex-ridge skip block entirely** (the PASS-2
   `convexRidge` detection). bracket/ahmed strips pass `saneCheck` → kept (recover
   +45 / +184); cow-ear strips fail → demoted → cow stays clean. No angle, no
   per-geometry tuning — it adapts to the geometry by construction.

This validates the **actual** built faces (no prospective-geometry duplication, so
no drift), which is why it generalizes to any degenerate-strip cause, not just
convex ridges.

## Hunks (anchors are code, not line numbers — the file shifted after the S0 patch)

### H1 — declare the demote state (next to the skip counters, after `nConvex`)

```cpp
    labelHashSet failedCols;                 // columns whose built faces fail saneCheck
    List<label> nDemoted(specs.size(), 0);   // construction-validated demotes (I1.2 rung-1)
```

### H2 — `saneCheck` collects the failing column (inside the lambda, where it emits INSANE)

After the existing `++nInsane; Info<< "SPLITLAYERS|INSANE|..." ...;` block, add:

```cpp
            failedCols.insert(ci2);          // drive a demote, not just a report
```

(Keep the `nInsane >= 20` early-return only for the *Info* spam cap; still record
the column even past 20 so the demote set is complete — move the cap to guard only
the `Info<<`, not the `failedCols.insert`.)

### H3 — remove the geometric convex-ridge skip (PASS 2)

Delete the whole `// --- convex-ridge opposed-column skip (I0.1 ...) ---` block
(the `convexRidge` detection that sets `col.skip="convexRidge"`). Leave the
`-convexRidgeAngle` option parsing in place (harmless, now unused) and keep
emitting `skippedConvexRidge|<nConvex>` (now always 0) for report back-compat; add
`|demoted|<nDemoted[si]>` to the patch report line.

### H4 — the rebuild loop (the delicate part)

Wrap the construction. Sketch:

```cpp
    const label maxBuild = 2;                 // one validate pass + one rebuild
    for (label buildAttempt = 0; buildAttempt < maxBuild; ++buildAttempt)
    {
        nInsane = 0;
        failedCols.clear();
        // reset shared state the construction mutates:
        forAllIters(railsOf, it) { it().rings.clear(); }   // PASS 3 refills these

        polyTopoChange meshMod(mesh);
        // ... PASS 3 ring creation ...
        // ... PASS 4 topology ...
        // ... lateral / absorption passes ...
        // ... (optimize is sweeps=0 by default; fine inside) ...

        if (failedCols.empty())
        {
            // sane → apply + write below (break out and proceed to write)
            break;
        }
        if (buildAttempt + 1 >= maxBuild)
        {
            break;   // fall through to the I4.1 self-check, which will refuse
        }
        // demote the offenders and rebuild once
        for (const label ci : failedCols)
        {
            if (!columns[ci].skip)
            {
                columns[ci].skip = true;
                columns[ci].reason = "demoted";
                ++nDemoted[columns[ci].si];
            }
        }
        // (meshMod is block-scoped → discarded; loop re-enters and rebuilds)
    }
    // I4.1 self-check (nInsane>0 → split-FAILED) then setInstance + write, as today.
```

**Care points (verify on WSL):**
- `meshMod` and the PASS-3/4 locals (`segRings`, `edgeClaim`, `newPtPos`,
  `stackOf`, `cellCol`, `nEdgeConflict`) must be **declared inside** the loop so
  each iteration is fresh. Move their declarations in.
- `railsOf[*].rings` is the one piece of shared state that persists across
  iterations — clear it at the top of each (shown). Confirm nothing else
  (e.g. `claimed`, the PASS-2 records) is mutated by PASS 3+; PASS 2 ran once,
  before the loop, so its skips persist (correct — demotes accumulate on top).
- The write (`mesh.setInstance/mesh.write`) and the final `SPLITLAYERS|result`
  stay **after** the loop (unchanged), guarded by the I4.1 self-check.

## DoD / what to measure (WSL, edge / sweeps 0, --checkmesh)

- **bracket**: Mesh OK, columns ≥ 1607 (recovers the +45), `demoted` small.
- **ahmed**: Mesh OK, columns ≥ 1891 (recovers the +184), `demoted` ~0.
- **cow**: Mesh OK at ≥ 98% coverage, `demoted` ≈ the old `skippedConvexRidge`
  (~240), `skippedConvexRidge`=0. The ears are demoted because their strips fail
  `saneCheck`, not because of an angle.
- lbend/cavity smoke: unchanged.
- Capture **bracket_corner** + **cow_ear** into the I4.2 matrix (cow_ear must fail
  on the pre-this-patch binary — it did, 30+40 faces — and pass after).

## Relationship to the rest of I1.x

This is **I1.2 rung-1** (demote→rebuild) + the **I4.1** orientation/quality
backstop, and it **retires the geometric `skippedConvexRidge`** (the I0.1
placeholder). It does NOT touch the conflict skips — the S0 data showed those are
scale-only (manifold/ahmed → S1 scale consensus) and joints/corner-bisector
(bracket/cow → S2 keying). S1 and S2 (the shared-rail consensus,
`2026-06-16-shared-rail-consensus.md`) are still the next levers for
`skippedConflict`. Later I1.2 rungs add the intermediate fallbacks (n−1, n−2,
single compressed layer) before demoting all the way to skip.
