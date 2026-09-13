# segqa — automated, label-free QA for published segmentations

**Repo:** https://github.com/Wadoekeani/vc-segqa  ·  Python + a browser viewer, no credentials, no bulk download

> [中文版](SUBMISSION.zh-TW.md)

## The problem

Segmentation quality is checked by eye, one segment at a time. Automating it
runs into a bootstrapping problem: validating a quality metric needs segments
already labelled good or bad, and no such labels are published. Every metric
you invent is unfalsifiable, so nobody builds one.

The way around it is to find questions that are **objectively answerable from
the published data alone**, where nobody has to agree with you about the
answer. There are two, and they turn out to catch different things.

## Check 1 — is the surface even on the scroll?

The published volumes are **masked**: voxels outside the scroll body are exactly
`0`. So "does this vertex land on the object at all" needs no labels, and a
surface that leaves the scroll is wrong by inspection.

Every published segment with a masked volume to check against — **299 segments,
11 scrolls, 147 M vertices**:

| scroll | segments | outside the scroll |
|---|---|---|
| PHerc0800 | 6 | **0.00 %** |
| PHerc0814 | 19 | 0.05 % |
| PHerc0172 | 53 | 0.14 % |
| PHerc0139 | 38 | 0.25 % |
| PHerc0343P | 8 | 3.26 % |
| PHerc0009B | 18 | 5.39 % |
| PHerc0500P2 | 38 | 5.45 % |
| PHercParis4 | 81 | 7.39 % |
| PHerc0841 | 3 | 11.82 % |
| PHerc1667 | 20 | **19.04 %** (17 of 20 segments over 5 %) |
| PHerc1447 | 15 | **33.40 %** |

Sharply bimodal. That four scrolls score at or near `0.00 %` is the argument
that this measures something real rather than boundary noise.

On Scroll 1 the rate rises monotonically with winding number — 0.00–7 % for
w010–w100, 16–39 % for w120–w129. The outermost segment leaves the object over
most of the upper scroll: 76 % of its vertices above z≈3700 sit in empty space.
Drawn against the CT, an entire arc of the trace runs through void:

![escape](outer-winding-escape-z3700.png)

## Check 2 — is it on the *right* sheet?

Check 1 is blind to a tracer that slips onto the neighbouring sheet: the
surface is still inside the scroll. So is distortion, because the
parameterisation stays locally isometric across a switch. The signal does not
exist inside a single segment.

It exists *between* them. Segment names carry official winding numbers, and two
segments one winding apart must sit one papyrus thickness apart everywhere they
overlap. That expectation is free, exactly like the mask.

Over the 92 adjacent pairs that exist directly:

| scroll | pairs | flagged |
|---|---|---|
| PHerc0172 | 43 | 0 |
| PHerc0139 | 35 | **2** |
| PHerc1667 | 14 | 0 |

**PHerc0139 w46 has jumped onto w45's sheet:**

```
w44/45   1.28   normal
w45/46   0.25   <-- coincident
w46/47   2.09   <-- double gap
w47/48   1.11   normal
```

Self-consistent, and it names the culprit: had w45 been the one that moved,
w44/45 would read double, and it does not. Rendered, the two surfaces
interpenetrate across their whole extent:

![sheet switch](viewer-sheetswitch-PHerc0139-w045-w046.png)

### Checks 1 and 2 are independent, and that is the point

PHerc1667 is the **worst** scroll by mask escape (19 %) yet is perfectly clean
here, with the tightest ratio range of the three — its problem is drifting off
the object at the outer edge, not confusing sheets. PHerc0139 is among the
**cleanest** by mask (0.25 %) and is the only one carrying a sheet switch.

A scroll can look flawless under one check while the real defect is visible
only to the other. Neither is sufficient alone.

## Check 3 — do two volumes agree about where the sheet is?

Scroll 1 now exists at four resolutions, and each fine volume publishes a
`transform.json` mapping its voxels onto the 7.91 µm volume from 2023. So most
segments have been meshed twice, and a shared frame is published too. A third
expectation follows, as free as the first two:

> two meshes of the same segment, cut on two different volumes, must describe
> the same sheet.

`xres.py` measures point-to-surface distance — nearest vertex, then the distance
to a plane fitted through its 3×3 grid neighbourhood. Point-to-*vertex* would
not do: the grids are ~158 µm apart, so querying a surface with points known to
lie on it already returns a median of 79 µm. Point-to-surface drops that floor
to 5.4 µm.

Across the 55 Scroll 1 segments meshed on both the 45.5 µm and the 7.91 µm
volume:

| pair | segments | disagreement p50 | = voxels of the coarser volume |
|---|---|---|---|
| 45.5 µm vs 7.91 µm | 55 | **46.6 µm** | 1.02 |
| 2.4 µm vs 7.91 µm | 9 | 3.7 µm | 0.47 |
| 1.129 µm vs 2.4 µm | 1 | 0.34 µm | 0.14 |

The lower two rows are the control arm, and they are what make the top row
readable: same comparison, same transforms, same code, twelve and a hundred
times tighter. The 46.6 µm is therefore not the method and not the transform
format — a surface traced at 45.5 µm is uncertain by about one voxel of that
volume.

![cross-resolution agreement](xres-agreement.png)

Two structures inside it. **Nine `5753_*` segments sit at 82–90 µm while the
other 46 sit at 33–64 µm** — no overlap, a whole batch flagged without anyone
labelling a vertex. And **disagreement rises outward**, p50 against winding
number at r = 0.86, from ~35 µm at w010–w045 to ~62 µm at w116–w129;
three winding ranges were segmented twice and the two batches agree to 2.3–4.6
µm, so the level reproduces. That matches the mask check's own outward trend
without being circular — one asks whether a vertex is on the object, the other
whether two volumes put the same sheet in the same place.

The control arm samples all three groups, not only the one that scores best at
45.5 µm. The two `5753_*` segments come back at 2.9 and 3.3 µm — *tighter* than
the 2023 batch — and `w128-129`, the worst of the winding trend at 45.5 µm,
returns 3.8 µm, the same as the innermost range. Both structures in the top row
therefore live on the 45.5 µm side: neither that batch nor the outer windings
are intrinsically hard to place, only hard to place *at 45.5 µm*. (w128-129 does
carry a longer tail even here — p90 29 µm against 11–14 elsewhere.)

What it does not resolve: the 45.5 µm side carries two candidate causes this
data cannot separate — the coarse volume's resolution limit, and the quality of
its own registration (its `transform.json` reproduces its landmarks to ≈356 µm,
where the 2.4 µm one manages ≈12 µm). Pure misregistration would appear as a
local shift and does not: inside 4 mm blocks the mean signed distance is ~28 µm
against a ~77 µm within-block spread, so most of the disagreement is shape.

No "fraction beyond half a sheet" is reported, though it would read better.
Every sheet-spacing figure available here is itself a nearest-point distance on
a ~900 µm grid with a discretisation floor of unknown size; dividing by it would
launder that into a clean-looking percentage. The coarse voxel is a denominator
that is actually known.

## Scroll 1, and one hypothesis eliminated

Scroll 1 names its segments by winding *range* (`w010-027`), so it had no
adjacent pairs. `unroll.py` recovers a per-vertex winding number from the angle
about the scroll axis, using the range in the name as calibration.

That calibration is itself the validation: **all 58 segments measure an exact
integer number of turns** (4.0, 3.0, 2.0 …), which a mis-fitted axis or a
mis-aligned unwrap could not produce. It also showed the naming is inclusive —
`w116-117` is two turns, not one.

With windings recovered, the same comparison runs: **230 adjacent pairs across
both batches, zero flagged**, ratios 0.68–1.32.

This settles a live question about the missing Scroll 1 title. The title sits
in the innermost windings, and all three geometric checks pass there: the core
*is* traced, its segments are 0.00–1.74 % outside the mask (thirteen of
twenty-one are exactly 0.00 %), and the winding topology is sound at w010–w031.

Check 3 agrees from its own direction: the innermost windings have the *lowest*
cross-resolution disagreement of any range (34–44 µm at w010–w045, against
60–64 µm at the outer edge), so the coarse and fine volumes place the core's
sheets in the same spot.

**The surfaces in the title region are geometrically healthy, so "no ink was
detected there" is not explained by a bad segmentation.** That points effort at
the ink-detection side rather than the tracing side. Ruling a hypothesis out is
less satisfying than an answer, but it says where the effort should not go.

### The core, asked directly

The core has no meshes to check, so `core_structure.py` asks the 2.4 µm image
instead, at twelve heights, with the traced band just outside w010 as the
control on the same slice. Two quantities: how much of the area has structure
(coverage), and whether that structure lines up into sheets (orientation
coherence over 96 µm, which uses no axis).

The sheets inside w010 are as well ordered as the traced ones — paired
coherence difference −0.006, unchanged under three thresholds. In the lower half
of the scan they are as plentiful too (40.6 % vs 38.7 %); in the upper half a
third less (25.8 % vs 42.3 %), and higher slice numbers are the top of the
scroll — the quantitative form of the Title Prize page's *"the top rows are
physically missing"*. A longitudinal cut shows the core's sheets continuous but
strongly inclined, one side of the core empty; cross-sections of the same
region read as fragments, which was an oblique-cut artefact, withdrawn.

So in the lower half of the scanned range the image does not explain why
tracing stopped at w010. In the upper half the papyrus is genuinely sparser.

## Controls

**Resolution.** Scrolls are sampled at different effective resolutions because
the pyramids stop at level 5. Rather than assume that is harmless, it was
measured: the same segments across an 8× range of voxel size move by 1.29×,
biased the way you would expect — a downsampled mask is blurred outward and so
more forgiving — with the ranking unchanged at every level. Against a 0 %–33 %
spread between scrolls it changes no conclusion. Absolute rates are ±30 %;
rankings are sound.

**Physical sanity.** The inter-sheet spacing is never supplied as a constant,
it falls out: 147 / 151 / 164 µm on three scrolls whose voxel sizes differ by
4×, computed independently. Scroll 1 reads higher at ~390 µm; it is a different
scroll and coarsely sampled, and the number is reported rather than explained.

**A negative result, kept.** The first approach — grid distortion on the tifxyz
quads — does not work, structurally rather than for want of tuning: across a
sheet switch the parameterisation stays locally isometric, so distortion is
blind to it by construction. On real data it fires only along hole boundaries,
which the validity mask already gives for free. Documented so the next person
does not spend a day rediscovering it.

**A second dead end, also kept.** Walking one turn along a grid row inside a
single segment measures the slant of the row, not the sheet gap: on Scroll 1's
outer windings the landing point is 2194 µm away in z against a 471 µm radial
step. Measuring only the radial component brings the median to a plausible
212 µm, but the within-segment spread is too wide to localise anything. Not
used; superseded by the unrolling above.

## Viewer

![viewer](viewer-screenshot.png)

A percentage says a segment is 30 % wrong. It does not say *where*, or in what
shape. The viewer draws traced surfaces in 3D with escaped regions in red — in
the shot above the entire upper scroll is red at a glance — and lists the
adjacent-winding pairs with their gap ratios. Click a flagged pair and it loads
exactly those two windings, coloured by identity, so a sheet switch is visible
as interpenetration rather than as a number.

All 11 scrolls and 299 segments. Meshes stream straight from the open-data
bucket (it is CORS-open); only per-vertex flag files are served locally,
subsampled to at most 128×128 per segment, so the whole viewer is 2.9 MB rather
than the 106 MB the 2.4 µm meshes would otherwise need.

## Notes that may be useful regardless

- `tifxyz` is three float32 TIFFs. 376 of 380 planes are classic, uncompressed,
  single strip — two lines to read. The other four are BigTIFF, tiled, LZW with
  the floating-point predictor, and a classic-only reader fails on them with a
  struct error that names nothing. Reading those four moved PHerc1447 from
  52 % to 33 %.
- The volumes are uncompressed, so one z-plane inside a chunk is a contiguous
  16 KB block. A Range request pulls a full-resolution cross-section for
  ~1.3 MB instead of 162 MB.
- A missing chunk is not an error: zarr's `fill_value` is `0` here, which means
  outside the scroll. 404 is an answer.
- Meshes are registered to a specific volume (`-on-<volume-id>-<res>um.tifxyz`).
  Sampling against a different volume of the same scroll gives silent garbage.
- Scroll 1's axis wanders ~400 voxels over its height; a single global centre
  makes radius meaningless near the core.
- Winding numbers run inside-out — w001 is innermost — and ranges in names are
  inclusive.

## Limitations

- Check 2 detects whole-winding jumps over a substantial area. A local switch
  affecting a small patch is diluted in the median and would be missed.
- Check 2 needs segments that overlap; 8 pairs in total were skipped for
  insufficient overlap, and those are reported rather than silently dropped.
- Check 1 needs a masked volume, which excludes 5 scrolls that have segments.
- "Outside the mask" is not always "the tracer was wrong": the other reading is
  that the mask excludes frayed outer layers that really are papyrus. Both
  matter, only the first is a bug, and this does not separate them in general.
  The rendered cross-sections settle it case by case.
- 12 of 311 published segments are not scanned: they carry no mesh registered
  to the scroll's masked volume.
