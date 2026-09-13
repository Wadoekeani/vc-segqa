# segqa — label-free quality checks for Herculaneum scroll segmentation

Automated QA for published `tifxyz` surface meshes. Reads everything straight
from the Vesuvius Challenge open-data bucket over plain HTTP — no credentials,
no bulk download, no local copy of a scroll.

> [中文版](README.zh-TW.md) · prize write-up: [SUBMISSION.md](SUBMISSION.md)

## The idea

Segmentation QA has a bootstrapping problem: to validate a quality metric you
need segments already labelled good or bad, and no such labels are published.
Every metric you invent is unfalsifiable.

There is one exception. The published volumes are **masked** — voxels outside
the scroll body are exactly `0`. So

> *does this vertex land on the object at all?*

is objectively answerable with no labels, for every vertex of every published
mesh. A mesh that leaves the scroll is wrong, and nobody has to agree with you
about it.

That is what this tool measures.

## What it found

Scanning Scroll 1 (`PHercParis4`, 81 segments, 4.1 M vertices):

| windings | vertices outside the scroll |
|---|---|
| w010–w100 | 0.00 – 7 % |
| w101–w115 | 10 – 14 % |
| w120–w129 | 16 – **39 %** |

Monotone in winding number, and exactly `0.00 %` on the inner windings — the
metric is silent where it should be silent. The outermost segment
(`w128-129`) leaves the object over most of the upper scroll: 76 % of its
vertices above z≈3700 sit in empty space.

![outer winding escape](outer-winding-escape-z3700.png)

Left: raw CT. Right: the traced surface. The left and top of the trace hug the
scroll edge correctly; the right-hand arc runs through pure void.

Two more results from the same data:

- **A real coverage gap.** At z≈300 (the bottom ~18 mm) the papyrus is intact
  and clearly layered, but only 8 traced rows cross it, against 86 at
  mid-height. Not damage — nobody has traced it.
  See `sections-z300-2000-3850.png`.
- **Version regression.** Each winding range was segmented twice (2026-06-23
  and 2026-07-01). The re-run improved the outer windings (w128-129: −8.9 pp)
  and regressed the inner ones (w046–w094: up to +2.6 pp, from a perfect
  0.00 %). It traded inner accuracy for outer reach.

## Corpus scan

299 segments across 11 scrolls, 147 M vertices:

| scroll | segments | vertices | outside | over 5% |
|---|---|---|---|---|
| PHerc0800 | 6 | 33 K | **0.00 %** | 0 |
| PHerc0814 | 19 | 1.4 M | 0.05 % | 0 |
| PHerc0172 | 53 | 28.6 M | 0.14 % | 0 |
| PHerc0139 | 38 | 4.4 M | 0.25 % | 0 |
| PHerc0343P | 8 | 112 K | 3.26 % | 1 |
| PHerc0009B | 18 | 1.2 M | 5.39 % | 10 |
| PHerc0500P2 | 38 | 531 K | 5.45 % | 7 |
| PHercParis4 | 81 | 4.1 M | 7.39 % | 27 |
| PHerc0841 | 3 | 136 K | 11.82 % | 1 |
| PHerc1667 | 20 | 106 M | 19.04 % | 17 |
| PHerc1447 | 15 | 255 K | 33.40 % | 10 |

The distribution is sharply bimodal: four scrolls are effectively zero
(0.00–0.25 %), two are badly off. That four scrolls can score `0.00 %` is the
argument that this measures something real rather than boundary noise.

**Resolution bias is measured, not assumed.** Pyramids stop at level 5, so the
effective sampling resolution varies with each volume's base voxel size — 76.8 µm
for PHerc1667 against 276–300 µm for most others. To find out whether that
matters, `sensitivity.py` re-runs the same segments across an 8× range of voxel
sizes:

| segment | 91 µm | 182 µm | 364 µm | 728 µm |
|---|---|---|---|---|
| w069-072 | 0.50 % | 0.04 % | 0.01 % | 0.01 % |
| 20230702185753 | 3.04 % | 1.41 % | 1.31 % | 1.38 % |
| w085-088 | 3.74 % | 3.31 % | 3.30 % | 3.38 % |
| w120-121 | 20.11 % | 18.86 % | 18.36 % | 17.75 % |
| w128-129 | 38.72 % | 32.33 % | 30.00 % | 29.28 % |

The mean moves 11.17 % → 8.64 % across that whole range: a **1.29× spread for
8× in voxel size**, biased the way you would expect — a downsampled mask is
blurred outward and so more forgiving. The ranking is unchanged at every level.

The corpus scan spans a 4× resolution range, narrower than what was tested
here, so the cross-scroll bias is bounded well below 1.29×. Against a spread of
0.00 % to 33 % between scrolls — two orders of magnitude — it does not change
any conclusion. Treat absolute values as ±30 % and rankings as sound.

## Sheet switches: a second, independent check

The mask check cannot see a tracer that slips onto the neighbouring sheet — the
surface is still inside the scroll. Neither can distortion, because the
parameterisation stays locally isometric across a switch. The signal does not
exist inside a single segment.

It exists *between* segments. Names carry official winding numbers, and two
segments one winding apart must sit one papyrus thickness apart everywhere they
overlap. That expectation is free, exactly like the mask:

- gap collapses to ~0 → both traced the same sheet; one of them jumped
- gap roughly doubles → a winding was skipped
- which segment moved is identified by checking its *other* neighbour

`sheetswitch.py` over the 92 adjacent-winding pairs that exist:

| scroll | pairs | flagged | ratio range |
|---|---|---|---|
| PHerc0172 | 43 | 0 | 0.64 – 1.46 |
| PHerc0139 | 35 | **2** | 0.25 – 2.09 |
| PHerc1667 | 14 | 0 | 0.85 – 1.21 |

PHerc0139 w46 has jumped onto w45's sheet:

```
w44/45   1.28   normal
w45/46   0.25   <-- coincident
w46/47   2.09   <-- double gap
w47/48   1.11   normal
```

Self-consistent, and it identifies the culprit: had w45 been the one that
moved, w44/45 would read double, and it does not.

**The two checks are independent, and that is the point.** PHerc1667 is the
worst scroll by mask escape (19 %) yet is perfectly clean here, with the
tightest ratio range of the three — its problem is drifting off the object at
the outer edge, not confusing sheets. PHerc0139 is among the cleanest by mask
(0.25 %) and is the only one carrying a sheet switch. A scroll can look
flawless under one check while the real defect is only visible to the other.

The measured spacing lands at 147 / 151 / 164 µm on the three scrolls —
independently computed, from volumes whose voxel sizes differ by 4×, with no
physical constant supplied anywhere.

Two pairs were skipped for insufficient overlap (reported, not silently
dropped).

### Scroll 1, unrolled

Scroll 1 names its segments by winding *range* (`w010-027`), so it had no
adjacent pairs. `unroll.py` recovers a per-vertex winding number from the angle
about the scroll axis, calibrated against the range in the name — and that
calibration is the validation: **all 58 segments measure an exact integer
number of turns** (4.0, 3.0, 2.0 …), which a mis-fitted axis could not produce.
It also showed the naming is inclusive: `w116-117` is two turns, not one.

With windings recovered, `scroll1_switch.py` runs the same comparison: **230
adjacent pairs across both batches, zero flagged**, ratios 0.68–1.32.

That settles a question about the missing Scroll 1 title, which sits in the
innermost windings. All three geometric checks pass there — the core *is*
traced, its segments are 0.00–1.74 % outside the mask (thirteen of twenty-one
exactly 0.00 %), and the winding topology is sound at w010–w031. The surfaces
are geometrically healthy, so "no ink was detected there" is not explained by a
bad segmentation.

### The core, asked directly

Everything above measures meshes. The core has none — the innermost traced sheet
is w010, at 0.5–1.6 mm from the axis depending on height — so the only thing
left to ask is the image: is there anything in there a tracer could have
followed? `core_structure.py` asks the 2.4 µm volume, ~2 MB per slice over HTTP
Range, at twelve heights (36–177 mm along the axis), with the traced band just
outside w010 as the control on the same slice:

- **coverage** — fraction of the area whose gradient energy clears a global
  threshold: how much of it has structure at all
- **coherence** — |mean unit double-angle gradient| over a 96 µm window: do
  those structures line up into sheets

| | core (inside w010) | traced band (w010 outward) |
|---|---|---|
| coherence, median paired difference | −0.006 (range −0.087 … +0.046, n = 10) | — |
| coverage, lower half of the scan | 40.6 % | 38.7 % |
| coverage, upper half of the scan | **25.8 %** | 42.3 % |

The sheets inside w010 line up exactly as well as the ones that were traced —
the coherence difference stays at zero under three different thresholds
(−0.006 / −0.006 / −0.007 at the 40th / 60th / 80th percentile). In the lower
half there is also just as much papyrus. In the upper half there is a third
less — and since the official FAQ puts higher slice numbers at the top of the
scroll (all three `transform.json` files carry a positive z→z coefficient, so the
2026 volumes keep that orientation), this is the quantitative form of the Title
Prize page's *"the top rows are physically missing"*.

![core structure](core-structure.png)

A longitudinal cut through the axis, `core-longitudinal.png`, shows why the
core looks the way it does in cross-section: its sheets are continuous over
2 mm of height but strongly inclined, crossing the fitted axis, and one side of
the core is empty. Cut obliquely, an inclined sheet looks like fragments; an
earlier reading of those cross-sections as "broken papyrus" was wrong and is
withdrawn.

For the missing title — which by convention sits in the core, and in PHerc. 172
was found exactly there in 2025 — this means: in the lower half of the scanned
range the core holds as much papyrus, as well ordered, as the part that has been
traced. Nothing in the image explains why tracing stopped at w010 there. The
upper half is genuinely emptier.

Limits: two-dimensional slices; the band averages depend on the fitted axis
(coherence itself does not); one longitudinal cut.

## Cross-resolution agreement: a third check

Scroll 1 now exists at four resolutions, and most segments have been meshed on
more than one of the volumes. Each fine volume ships a `transform.json` mapping
its voxels onto the 7.91 µm volume from 2023, so a shared frame is published
too. That gives a third free expectation, in the same spirit as the first two:

> *two meshes of the same segment, cut on two different volumes, must describe
> the same sheet.*

No labels, no opinion about what a good mesh looks like — just two independent
descriptions of one physical surface that have to agree.

`xres.py` measures point-to-**surface** distance: nearest vertex, then the
distance to a plane fitted through that vertex's 3×3 grid neighbourhood. Nearest
*vertex* distance does not work here, and the control in `xres.py` says why —
the grids are ~158 µm apart, so querying a surface with points known to lie on
it already returns a median of 79 µm. Point-to-surface brings that floor down to
5.4 µm, which is the measurement's own noise.

`xres_scan.py` over the 55 Scroll 1 segments that exist on both the 45.5 µm and
the 7.91 µm volume, with `xres_control.py` as the control arm:

| pair | segments | disagreement p50 | = voxels of the coarser volume |
|---|---|---|---|
| 45.5 µm vs 7.91 µm | 55 | **46.6 µm** | 1.02 |
| 2.4 µm vs 7.91 µm | 9 | 3.7 µm | 0.47 |
| 1.129 µm vs 2.4 µm | 1 | 0.34 µm | 0.14 |

The control arm is what makes the first row readable. Between the fine volumes
the same comparison, the same transforms and the same code return 3.7 µm and
0.34 µm — twelve and a hundred times tighter. So the 46.6 µm is not the
measurement, and not the transform file format. A papyrus sheet is only a few
45.5 µm voxels thick, and a surface traced at that resolution turns out to be
uncertain by roughly one of them.

Two structures inside that number:

**The 5753 family is twice as bad as everything else, with no overlap.**

| group | segments | p50 range | median, in 45.5 µm voxels |
|---|---|---|---|
| `w`-named | 33 | 33 – 64 µm | 0.97 |
| everything else (11 from 2023, 2 re-meshes) | 13 | 39 – 61 µm | 0.96 |
| `5753_*` family | 9 | **82 – 90 µm** | **1.90** |

Nine segments, every one of them worse than all 46 others. That is the kind of
thing this check is for: it points at a specific batch without anyone having to
label a single vertex.

**Disagreement rises monotonically outward.** Over the `w`-named segments,
p50 correlates with winding number at r = 0.86:

| windings | p50 | p90 |
|---|---|---|
| w010–w045 | 34 – 44 µm | 86 – 151 µm |
| w046–w088 | 33 – 42 µm | 93 – 141 µm |
| w089–w115 | 46 – 52 µm | 160 – 184 µm |
| w116–w129 | 60 – 64 µm | 180 – 207 µm |

![cross-resolution agreement](xres-agreement.png)

Left: every `w`-named segment, coarse-vs-fine disagreement against winding
number, with the control arm as the two dotted floors. Right: the same p50 by
group — the `5753_*` family separates cleanly from everything else.

Three of the ranges were segmented twice (2026-06-23 and 2026-07-01); the two
batches agree to 2.3, 3.7 and 4.6 µm there, so the level is reproducible and not
an artefact of one run.

The control arm samples all three groups, not only the one that scores best at
45.5 µm. The two `5753_*` segments come back at 2.9 and 3.3 µm — *tighter* than
the 2023 batch — and `w128-129`, the worst of the winding trend at 45.5 µm,
returns 3.8 µm, the same as the innermost range. Both structures in the top row
therefore live on the 45.5 µm side: neither that batch nor the outer windings
are intrinsically hard to place, only hard to place *at 45.5 µm*. (w128-129 does
carry a longer tail even here — p90 29 µm against 11–14 elsewhere.)

The outward trend also matches the mask check, which rises outward on the same
data, and the agreement is not circular: one check asks whether a vertex is on the object, the
other asks whether two volumes place the same sheet in the same spot.

**What this does not resolve.** The 45.5 µm side carries two candidate causes
that this data cannot separate: the coarse volume's own resolution limit, and
the quality of *its* registration — the 45.5 µm `transform.json` reproduces its
own landmarks to 45 fixed-frame voxels (≈356 µm) where the 2.4 µm one manages
1.5 (≈12 µm). A pure misregistration would show up as a local shift, and it does
not: inside 4 mm blocks the mean signed distance runs ~28 µm against a
within-block spread of ~77 µm (median over the 11 segments compact enough to
fill three such blocks), so most of the disagreement is shape, not offset.
That argues against registration being the whole story without ruling it out.

**No "fraction beyond half a sheet" is reported**, though it would read well.
Every sheet-spacing number available here is itself a nearest-point distance on
a ~900 µm grid, carrying a discretisation floor of unknown size; dividing by it
would launder that uncertainty into a clean-looking percentage. The voxel of the
coarse volume is a denominator that is actually known.

## Viewer

![viewer](viewer-screenshot.png)

A Three.js view of the traced surfaces with the escaped regions in red, so the
*shape* of a failure is visible rather than just its percentage. In the shot
above the whole upper part of the scroll is red — the same defect that reads as
"76 % of vertices above z≈3700" in the table.

```bash
python3 export_flags.py scan_plan.json            # all scrolls -> viewer/data/
cd viewer && python3 -m http.server 8731
```

All 11 scrolls and 299 segments, picked from a dropdown. Meshes stream straight
from S3 (the bucket is CORS-open); only per-vertex flag files are served
locally — subsampled to at most 128×128 per segment, so the whole set is 2.2 MB
rather than the 106 MB the 2.4 µm meshes would otherwise need. Runtime state is
on `window.__segqa`.

The sidebar also lists the adjacent-winding pairs with their gap ratios. Click
one and it loads exactly those two windings, coloured by identity:

![sheet switch](viewer-sheetswitch-PHerc0139-w045-w046.png)

PHerc0139 w045 (green) and w046 (orange) should sit one papyrus thickness
apart. They interpenetrate across the entire surface — the ratio of 0.25 made
visible.

## Usage

Needs only `numpy` and `matplotlib`.

```bash
python3 survey.py                      # which scrolls have segments + a masked volume
python3 plan_scan.py PHercParis4 ...   # pick a volume and pyramid level per scroll
python3 scan.py scan_plan.json         # the scan -> scan_results.json
python3 report.py scan_results.json    # summary tables + figures
```

Downloaded meshes and volume chunks are cached under `~/.cache/segqa`
(override with `SEGQA_CACHE`). It is a cache, not a source — deleting it costs
a re-download, nothing else.

Cross-sections at any height and resolution, with the traced surfaces drawn on
top:

```bash
python3 section.py ~/.cache/segqa/PHercParis4 1 300 2000 3850 sections.png
```

Cross-resolution agreement needs the second mesh of each segment plus the
published transforms:

```bash
python3 xres_fetch.py 7.91um      # transforms + every 7.91 um mesh (~180 MB)
python3 xres_fetch.py 2.4um       # the control arm's other side (~2 GB; a few
python3 xres_fetch.py 1.129um     #   segments are enough - interrupt when bored)
python3 xres.py 20230702185753    # one segment at three resolutions, with the control
python3 xres_scan.py              # all of them -> xres_results.json
python3 xres_control.py           # control arm: the fine volumes against each other
python3 xres_plot.py              # -> xres-agreement.png
python3 core_structure.py         # the untraced core, asked directly -> core-structure.png
```

## Notes on the data

Things that cost time to discover, in case they save you some:

- **`tifxyz` is nearly trivial to read — until it isn't.** Of Scroll 1's 254
  published meshes, 240 are classic TIFF, uncompressed, single strip: two lines
  of code. The other 14 are BigTIFF, tiled, LZW with the floating-point
  predictor, and 8 of those are split across 2–12 tiles, which means the tile
  offsets live in an array the value field only points at. A single-tile reader
  fails on them with `buffer is smaller than requested size`, naming nothing.
  Resolution is no guide: 78 of the 81 2.4 µm meshes are plain single strips
  while 11 of the 55 coarser 7.91 µm ones are tiled (counted on each mesh's
  `x.tif`, read over HTTP Range — two requests per mesh, no download). `l0.py` handles all of it
  and still has no TIFF dependency.
- **Meshes are registered to a specific volume.** The directory name says
  which: `<seg>-on-<volume-id>-<resolution>um.tifxyz`. Sampling a mesh against
  a different volume of the same scroll gives silent garbage.
- **The volumes are uncompressed.** A chunk is exactly `128³` bytes in C
  order, so one z-plane inside it is a contiguous 16 KB block: a Range request
  fetches a full-resolution cross-section for ~1.3 MB instead of 162 MB.
- **A missing chunk means empty.** Zarr fills absent chunks with `fill_value`,
  which is `0` here — outside the scroll. So a 404 is an answer, not an error.
- **The scroll is bent.** Its axis wanders by ~400 voxels over the height, so a
  single global centre makes radius meaningless near the core. Fit per z-slab
  (`axis.py`).
- **Winding numbers are in the segment names** (`w010-027`, `-w062`) and run
  **inside-out**: w001 is the innermost. Verified — median radius rises
  monotonically across all 28 winding groups of Scroll 1.

## What does not work

`l0.py` measures grid distortion on the `tifxyz` quad mesh — edge anisotropy
and stretch. **It cannot detect sheet switches**, and the reason is structural,
not a matter of tuning: when a tracer walks onto the neighbouring sheet the
parameterisation stays locally isometric — cells keep their ~20 voxel spacing
and stay square. Distortion is blind to the defect by construction.

On real data it fires only along the ragged boundary of holes, which the
validity mask already gives you for free:

![L0 false positives](l0-false-positives.png)

The interior of every published segment is uniformly clean (p95 stretch 1.05 –
1.11). Kept here as a negative result so the next person does not spend a day
rediscovering it.

`wrapgap.py` is a second dead end, also kept. Walk one full turn along a grid
row and measure the distance to where you land, which should be one sheet
spacing. Measuring that in 3D is wrong: grid rows slant in z, badly on the
outer windings, where the landing point is 2194 µm away in height against a
471 µm radial step — so the 3D number measures the slant. Taking only the
radial component brings the median to a plausible 212 µm, but the spread within
a segment is too wide to localise anything. Superseded by `unroll.py`.

## Limitations

- **Sheet switching is detected at whole-winding scale.** A local switch over a
  small patch is diluted in the median and would be missed.
- **Check 2 needs overlapping segments**; 8 pairs in total were skipped for
  insufficient overlap, and those are reported rather than silently dropped.
- **"Outside the mask" is not always "the tracer was wrong."** The other
  reading is that the mask excludes frayed, detached outer layers that really
  are papyrus. Both matter, but only the first is a bug, and this tool does not
  separate them. The rendered cross-sections do: in the Scroll 1 case the trace
  runs through visibly empty space, which settles it there but not in general.
- **Needs a masked volume**, which excludes 5 scrolls that do have segments.
- **Absolute rates carry ~±30 % resolution uncertainty**; rankings are sound.

## Status

Working, not polished. Self-checks live in `l0.py`, `scan.py` and
`export_flags.py` — `python3 l0.py` runs its own with no data at all — but
there is no test suite beyond those.
