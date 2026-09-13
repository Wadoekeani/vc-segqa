# segqa — label-free quality checks for Herculaneum scroll segmentation

Automated QA for published `tifxyz` surface meshes. Reads everything straight
from the Vesuvius Challenge open-data bucket over plain HTTP — no credentials,
no bulk download, no local copy of a scroll.

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

295 segments across 11 scrolls, 147 M vertices:

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
| PHerc1447 | 11 | 162 K | 52.42 % | 10 |

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
0.00 % to 52 % between scrolls — two orders of magnitude — it does not change
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
dropped). Scroll 1 cannot use this check yet: its segments are named by winding
*ranges* (`w010-027`), so no adjacent single-winding pairs exist.

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

All 11 scrolls and 295 segments, picked from a dropdown. Meshes stream straight
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

## Notes on the data

Things that cost time to discover, in case they save you some:

- **`tifxyz` is trivial to read.** Three uncompressed single-strip float32
  TIFFs. The whole reader is 12 lines (`l0.py`); no TIFF library needed.
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

`wrapgap.py` is an unfinished second idea: walk one full turn along a grid row
and measure the 3D distance to where you land, which should be one sheet
spacing. On the one segment with enough angular coverage it recovers 294 µm —
the known papyrus spacing — without being told any physical constant. But
segments narrower than one full turn cannot be tested this way, so it is
unvalidated. Recorded, not relied on.

## Limitations

- **Sheet switching needs adjacent windings.** The mask check alone cannot see
  it; `sheetswitch.py` can, but only where segments carry single winding
  numbers. Scroll 1 names its segments by winding range, so it is not covered
  yet.
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
