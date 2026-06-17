# splitLayers: task 19 — OPPWALL instrument for the opposing-wall gap collision (2026-06-17)

**Track:** Phase 5 / boundary-layer engine — internal-flow (`manifold_gap`).
**Status:** instrument only (behaviour-neutral, emit-only). **No construction change.**
This is the diagnostic-first step the handoff §3 demands: replace the opposing-wall
*hypothesis* with *geometry* before designing any fix.

## Why an instrument, not a fix

`manifold_gap` fails checkMesh with **17790 open cells** (current tip: convex skip
retired to 179 + gap-aware variable-`n`). The three prior reads were wrong or
incomplete:

1. ❌ variable-`n` termination (deeper-defer) — barely moved it (18282→17790).
2. ❌ convex-ridge opposed strips (the I0.1 cow class) — **disproven** by `ridgePts`:
   ahmed/bracket have 171/314 ridge points and **zero** open cells; the manifold has
   17790 open cells but only 1876 ridge points, not co-located. (So the premise of
   `2026-06-17-i13-opposed-strip-ridge-consensus.md` does NOT hold for the manifold.)
3. ◐ the gap-budget count-cap (G2/G3) — only **partially** relieved it; 17790 remain,
   so the failure is not purely a layer-count budget.

The surviving hypothesis is the **opposing-wall gap collision**: columns grown from
the two facing gap walls meet in the gap interior and their `seg0` strips
(`strips:seg0:cuts8:ocNNN:oseg0`) fail to close. The orchestrator cannot fix this
(it only demotes; demoting one wall opens the other → it exhausts). But we will NOT
write construction code on a hypothesis — we measure the geometry first.

## What the instrument adds

The existing `OPENCELL` attributor (after `changeMesh`) reports, per open cell, the
construct tags + the (first) owning column + its `nCol`. It does **not** report the
*coordinates* of the two facing columns' ring points — which is exactly what tells
us *how* the strips fail (overlap? offset? wrong winding?).

`OPPWALL` extends it to dump, for each opposed open cell (one whose faces name an
opposing column), the **3D coordinates of both columns' `seg0` ring points** on
their shared face, plus each rail's wall base point.

Everything needed is already in scope at the `OPENCELL` block:
- `newPtPos` (`Map<point>`, L1047/L1143) — ring point label → coordinate (points are
  not moved by `changeMesh`; `nOptSweeps == 0`).
- `columns` (rails / chain / `nCol`), `pts` (wall points), `claimed`, and the
  sharer-aware `ringsFor` / `assignFor` / `railKey` lambdas.

The only new state is `Map<label> faceOtherColOf` — origin face → the opposing
column `claimed[oth]` — populated next to the existing strip `recordFace` so the
attributor names **both** facing columns directly (no fragile tag-string parsing).

## The three hunks (all additive, `debugDump`-gated, emit-only)

1. **Declare** `Map<label> faceOtherColOf;` beside `faceTagOf`/`faceColOf` (~L258).
2. **Populate** it in the strip side-face pass, inside the existing `if (debugDump)`
   tag block, where `oc = claimed[oth]` is already computed (~L1768):
   `faceOtherColOf.insert(fi, oc);`.
3. **Emit** in the `OPENCELL` per-cell block (~L1903): a small `dumpSeg0(X, side)`
   lambda (guards skip / `railsOf.found` / `newPtPos.cfind` so it can never throw),
   and, for the matched `(colA, colB)` pair found on a face carrying both
   `faceColOf` and `faceOtherColOf`, emit:
   ```
   SPLITLAYERS|OPPWALL|cell|<celli>|openness|<x>|colA|<A>|nA|<n>|colB|<B>|nB|<n>
   SPLITLAYERS|OPPWALL|A|<A>|rail|<c>|wall|(x y z)|seg0|<k>:(x y z)<k>:(x y z)...
   SPLITLAYERS|OPPWALL|B|<B>|rail|<c>|wall|(x y z)|seg0|<k>:(x y z)...
   ```
   Bounded to the same first ≤16 open cells as `OPENCELL`.

Behaviour neutrality: all three are emit-only and gated on `debugDump`; no face,
point, cell, scale, or skip decision changes. ahmed/bracket/cow stay
**byte-identical** (the regression gate). Consistent with the existing
`OOB`/`badFace`/`UNHANDLED`/`OPENCELL` diagnostics — keep it on as a cheap guard.

## How the geometry forks the fix (what we look for in the output)

For each opposed open cell, compare `colA`'s and `colB`'s `seg0` ring coordinates and
the wall-to-wall gap:

- **Different positions, each its own ring set** → the shared gap-interior face is
  single-valued for one column only → fix = a cross-wall shared-seam / consensus on
  the gap face (the opposing-wall analogue of S2b's per-segment ring sharing).
- **Overlap / cross the gap midline** → stacks still too tall for the gap → residual
  budget; tighten the G3 margin / core band (or cap harder). Means the cap
  under-corrected, not a topology bug.
- **Coincident but opposite winding** → orientation bug in the opposed side-face
  owner/neighbour assignment when both sides are split columns from facing walls.

Each routes to a different, targeted fix — the reason to measure first.

## Build + run (David, WSL — `zefra-foam`, OF v2512)

NB the C++ app and the Python eval live in DIFFERENT repos — build from
`mesh-tools-gpl/`, run the eval from `Zefra/`.

```
# build (from the GPL repo — this is where splitLayers/Make lives)
cd /mnt/c/Users/traca/Desktop/fullCrossflow/mesh-tools-gpl
wmake applications/splitLayers

# eval (from the Zefra repo — the eval calls splitLayers off $PATH/$FOAM_USER_APPBIN)
cd /mnt/c/Users/traca/Desktop/fullCrossflow/Zefra
# 0. baseline: confirm the 17790 starting point + externals clean
python3 tools/2026-06-12-spike-split-eval.py all --placement edge --sweeps 0 --checkmesh
# 1. the instrument run
python3 tools/2026-06-12-spike-split-eval.py manifold_gap --placement edge --sweeps 0 --debug --checkmesh
```

Send back the `SPLITLAYERS|OPPWALL|…` lines (plus the existing `OPENCELL` and `gap`
lines). Confirm ahmed/bracket/cow are byte-identical (neutrality check).

## DoD for the instrument (not the fix)

- `OPPWALL` lines emitted for the manifold's opposed open cells with real
  coordinates; ahmed/bracket/cow byte-identical.
- From the coordinates we can state the actual failure mode (offset / overlap /
  winding) and pick the fix family — then a separate design doc + approval before any
  construction C++.

---

## Update 2026-06-17 — OPPWALL result (hypothesis REFUTED) + GAPDECIDE step

Ran on the manifold (16 open cells sampled). Parsed coordinates (16/16):

- **NOT opposing-wall strip collision.** Every paired column is an **edge-adjacent
  neighbour on the same wall** (shares ≥2 wall vertices; normals only 0.2°–21° apart).
- **Shared seg0 strips coincide EXACTLY** (ring Δ = 0.00e+00). The §3 premise — "their
  seg0 strips fail to close" — is **false**; those strips are byte-identical from both
  columns and close. The imbalance is a *different* face.
- The columns are **single-cell** (chain=1), **8 layers crammed into one ~9.6e-5-thick
  cell at scale ~0.06**, sliver aspect ~33.
- Each gap cell carries **two `wall+rings` faces** — a wall on **both** sides: the gap
  is **one background cell thick**. Open cells come in pairs per column: a `walls=2`
  cell (openness ~0.58) + a `walls=1` partner (~0.8–0.99).

**Corrected mechanism.** It is a gap defect, but not via colliding strips. The gap is
one cell thick; a column claims that cell and, because `col.scale = avail/needed`
(L511), its `nCol` layers span the **entire** cell depth → the remainder cell is
**zero-thickness** and its far face is the opposing wall → it can't close. These
columns came out **`nCol=8`, uncapped** → the anti-parallel gap-cap (`wn·on < -0.5`,
L505) returned **false** for them. Why it returned false is the open question.

(NB owner cells with >1 patch face are skipped as corners at L374, so the two wall
faces appear at *construction* in the gap cell, not in the owner cell.)

### GAPDECIDE — confirmatory instrument (emit-only, this patch)

PASS-1 gap block now emits, for each column that keeps FULL layers while heavily
scaled (`nCol==nLayers && scale<0.5` = the open-cell danger zone), bounded to 60 +
a `gapDangerTotal`:

```
SPLITLAYERS|GAPDECIDE|col|..|chain|..|blockCol|..|blockFace|..|topInternal|..|
  oppFace|..|oppPatch|..|wn.on|..|avail|..|needed|..|cumLast|..|isGap|..|nCol|..|scale|..
```

It answers WHY `isGap` is false (no `oppFace` found ⇒ `wn.on` sentinel `2.0`; vs
`oppFace` found but `wn.on ≥ -0.5`), and confirms the zero-remainder arithmetic
(`avail ≈ scale·cumLast`). Externals never trigger (scale ~0.7–0.85) ⇒ `gapDangerTotal 0`
there = the neutrality check.

**Run:** `wmake` (from `mesh-tools-gpl`), then
`python3 tools/2026-06-12-spike-split-eval.py manifold_gap --placement edge --sweeps 0 --debug --checkmesh`
(from `Zefra`). The `GAPDECIDE` lines decide the fix:
- `oppFace = -1` (not found) ⇒ fix = make the discriminator detect the one-cell gap
  (the opposing wall is the column's own `topFace`/`iface`, a boundary face).
- `oppFace ≥ 0` but `wn.on ≥ -0.5` ⇒ the opposing-normal test is too strict / wrong
  sign for this geometry.
Either way the cap must also **leave a remainder margin** (it already uses
`avail*gapMargin` in the gap branch; the bug is reaching that branch).

---

## GAPDECIDE result + the fix (2026-06-17)

`GAPDECIDE` confirmed it: every manifold danger column has `oppFace ≥ 0` on the
layered patch but `wn.on = +0.09 .. +0.40` (never `< -0.5`), so `isGap=0` and
**531019 / 538941 columns (98.5%)** keep all 8 layers crammed into a sub-mm depth.
The gap is a thin, **oblique/curved** channel — the anti-parallel test is the wrong
tool.

**Why a flag, not a geometric auto-detector.** Running `GAPDECIDE` on the externals
settled the gating question:

| geom | `gapDangerTotal` | `wn.on` of danger cols | `scale` | clean? |
|---|---|---|---|---|
| ahmed | 0 | — | — | yes |
| cow | 0 | — | — | yes |
| bracket | **53** | 0.087 – 0.80 | 0.146 – 0.33 | **yes** |
| manifold | **531019** | 0.09 – 0.40 | 0.063 – 0.20 | no |

ahmed/cow have **zero** danger columns (safe automatically). But bracket's 53 (its
reentrant corner) are **per-column indistinguishable** from the manifold's — their
`wn.on` (0.087, 0.090, 0.099, 0.111…) and `scale` (0.146) sit inside the manifold's
ranges. No `wn.on` / `scale` / `chain` threshold separates 53 from 531019. The real
difference is **flow type** (bracket external, manifold internal) — known to the
planner, not to local geometry. So the cap is **gated on an internal-flow flag**;
externals never set it and are byte-identical *by construction*. (A purely geometric
detector that also handles external thin gaps — slats — is a later refinement; not on
the external-aero-first MVP path.)

**The fix (implemented, behind `-gapCap`).** In the PASS-1 gap block, after the
existing anti-parallel `isGap`:

```cpp
if (gapCap && cumHeights[si].last() > avail + SMALL) { isGap = true; }
```

When `-gapCap` is set, any column whose full stack can't fit the available inward
depth takes the existing cap branch (nFit full-size layers within `avail*gapMargin`,
leaving a remainder), regardless of opposing-wall angle. `-gapCap` off ⇒ the block is
a no-op ⇒ externals byte-identical. Manifold caps to ~1–2 layers (avail ≈ 1.3e-4 fits
~1 full y1) + a healthy remainder.

For the spike the eval passes it via `--gap-cap` (manifold only); production wiring is
a per-patch `gapAware` field the planner sets in the layer spec (additive, v0.2-style).

**Verify (David, WSL):**
```
cd /mnt/c/Users/traca/Desktop/fullCrossflow/mesh-tools-gpl && wmake applications/splitLayers
cd /mnt/c/Users/traca/Desktop/fullCrossflow/Zefra
# externals: NO flag -> must stay byte-identical / Mesh OK
python3 tools/2026-06-12-spike-split-eval.py all --placement edge --sweeps 0 --checkmesh
# manifold: WITH the flag -> expect openCells crash down, gap cols capped to 1-2
python3 tools/2026-06-12-spike-split-eval.py manifold_gap --placement edge --sweeps 0 --gap-cap --debug --checkmesh
```

**DoD.** ahmed/bracket/cow byte-identical + Mesh OK (flag off); manifold `openCells`
→ ~0 (or orchestrator-converges, not exhausts) with `cappedCols` ≈ the gap columns
and layers present where the gap allows. If a residual remains, the quality-retry
orchestrator (driver) cleans it; if it's large, OPPWALL/OPENCELL stay on to re-diagnose.

### Result (WSL, 2026-06-17) — gating + budget WORK; variable-n topology exposed

- **Externals byte-identical (flag off):** ahmed 254050/268477, bracket 466957/491990
  (recovered 131), cow 241414/271082, all `cappedCols 0`, Mesh OK. Gating guarantee
  holds — zero external regression.
- **Budget correct (manifold, `--gap-cap`):** `cappedCols 569255` (was 38236),
  `scaleLowered 0`, `meanScale 1.0` — the cramming is gone; every kept layer is
  full-size. Diagnosis + discriminator validated.
- **NEW blocker = variable-n termination** (the handoff's flagged risk). The wide
  spread of capped counts (many 1–2, 37543 skipped to 0) breaks the lateral/side-face
  topology: `badFaces 388`, `deferredTop 7425` (was 740), and `BADFACE|lateral` lines
  with **both** owner & neighbour dangling (`danglOwn 1 danglNei 1`, `internalFi 1`)
  → `changeMesh` crash, exit 1. Pre-cap the manifold completed (17790 open, no crash);
  capping trades the open cells for a variable-n topology crash.
- **Lead:** `deferredTop` 740→7425 suggests the deeper-neighbour defer (side-face
  pass) defers faces the deeper column never builds → they fall to the lateral pass,
  which keeps the (removed) original owner/neighbour → dangle. Extended `BADFACE` to
  emit `ownCol/ownN/neiCol/neiN`; `ownN!=neiN` ⇒ defer-not-picked-up (fix = build the
  deferred face / remap in the lateral pass); `ownN==neiN` ⇒ a different gap. Re-run
  `manifold_gap --gap-cap --debug` and read the nCol jumps before the topology fix.

### Variable-n cause confirmed + the smoothing fix (2026-06-17)

`BADFACE` (with `ownN/neiN`): **all** dangling lateral faces are `nCol` **8 / 2**
pairs — a full-depth column abutting a capped-to-2 neighbour. So the 388 failures are
purely a **layer-count gradient**: the `term`-clamp closes gentle steps (so the old
38k-column cap never crashed) but not a jump of 6. Not a random topology bug.

**Fix — gradual layer-count termination (`gapSmooth`, gapCap only).** A relaxation
pass between PASS 1 and PASS 2: iterate the layered-patch face adjacency
(`pp.faceFaces()`) and lower each column's `nCol` to at most `min(split-neighbour
nCol) + 1`, to a fixpoint. The 8/2 cliffs become 8->..->2 gradients every transition
can close. Monotone-decreasing => terminates; min-plus => order-independent fixpoint
(determinism for the I4.3 print). Only **split/split** jumps are smoothed; split/skip
(nCol 0) is left alone (the lateral pass handles it, no dangle). Nothing is pushed to
0 (min split neighbour >= 1 => floor 1). Emits `SPLITLAYERS|gapSmooth|maxStep|..|
lowered|..|sweeps|..`. Externals never set `-gapCap` => no-op => byte-identical.

**Verify:** `wmake`, then
`manifold_gap --placement edge --sweeps 0 --gap-cap --debug --checkmesh`. Expect
`badFaces 0` (no dangling), `deferredTop` back small, `openCells` -> ~0 (or the
quality-retry driver converges). If a residual `badFaces` remains at `ownN/neiN`
differing by 1, the deeper transition-face construction (not the gradient) needs work.
Then the shipping driver `tools/2026-06-16-splitlayers-quality-retry.py manifold_gap`
(needs a `-gapCap` passthrough — add if the plain run is close but not clean).

### Result + never-crash guard (2026-06-17)

`gapSmooth` fired (`lowered 4988, sweeps 3`) and cut the dangling faces **388 → 40**
(99.8% of the 18106 transition defers resolve). Externals byte-identical again
(ahmed 254050/268477, bracket 466957/491990, both Mesh OK; cow's fail is its
pre-existing 58 convex-ear cells). So the gradient was the dominant cause; the
remaining **40** are deeper "defer-not-picked-up" transition faces that escape the
side-face pass even at ±1 steps.

A *dangling* face is a hard `changeMesh` FATAL, not a demotable cell — so the residual
can't go straight to the orchestrator. Added a **lateral-pass remainder-remap**
(`latRemap`, mirrors the top-face remap): any lateral face left with a removed chain
cell as owner/neighbour is reassigned to that column's remainder → a valid (if poor)
face instead of `(owner, -1, -1)`. The residual cells become open/skew → the
quality-retry orchestrator demotes the ~40 columns (well within its budget; it cleared
23–55 on the cow). No-op for unsplit owners ⇒ externals byte-identical (their lateral
faces have no removed-cell owners — that's why `badFaces 0` there).

Expected next run (manifold `--gap-cap`): `badFaces 0`, `result split-ok` (no crash);
then `2026-06-16-splitlayers-quality-retry.py manifold_gap` (with a `-gapCap`
passthrough) converges by demoting the residual. A clean *build* of those 40
transition faces (no demote) is a future refinement; the manifold ships on extrusion
meanwhile, so this is not on the external-aero MVP path.
