# segqa — automated, label-free QA for published segmentations

**Repo:** <fill in>  ·  **Live viewer:** run `viewer/` locally, see README

## The problem

Segmentation QA is checked by eye, one segment at a time. Building a tool to
automate it runs into a bootstrapping problem: validating a quality metric
needs segments already labelled good or bad, and no such labels are published.
Every metric you invent is unfalsifiable, so nobody builds one.

## The way around it

The published volumes are **masked** — voxels outside the scroll body are
exactly `0`. So this question:

> does this vertex land on the object at all?

is objectively answerable, for every vertex of every published mesh, with no
labels and no room for disagreement. A surface that leaves the scroll is wrong.

That single check turns out to carry a lot of signal.

## What the scan found

Every published segment with a masked volume to check against — **295 segments
across 11 scrolls, 147 M vertices**:

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
| PHerc1447 | 11 | **52.42 %** |

Sharply bimodal. Four scrolls are effectively zero; two are badly off. That
four scrolls can score `0.00 %` is the argument that this measures something
real rather than boundary noise.

### Three concrete results

**1. Outer windings drift off the object.** On Scroll 1 the escape rate rises
monotonically with winding number: 0.00–7 % for w010–w100, 16–39 % for
w120–w129. The outermost segment leaves the scroll over most of the upper
half — 76 % of its vertices above z≈3700 sit in empty space. Rendered against
the CT, an entire arc of the trace runs through void:

![escape](outer-winding-escape-z3700.png)

**2. An untraced region that is not damage.** At the bottom of Scroll 1
(z≈300) the papyrus is intact and clearly layered, but only 8 traced rows
cross it against 86 at mid-height. Nobody has segmented it.

**3. A measurable regression between batches.** Each Scroll 1 winding range
was segmented twice. The July re-run improved the outer windings (w128-129:
−8.9 pp) and regressed the inner ones (w046–w094: up to +2.6 pp, from a
perfect 0.00 %) — it traded inner accuracy for outer reach. This is the kind
of comparison the tool is meant to make routine.

## Controls

**Resolution.** Scrolls are sampled at different effective resolutions because
the pyramids stop at level 5. Rather than assume that is harmless, it was
measured: the same segments across an 8× range of voxel size move by 1.29×,
biased the way you would expect (a downsampled mask is blurred outward and more
forgiving), with the ranking unchanged at every level. Against a 0 %–52 %
spread between scrolls, that bias changes nothing.

**A negative result.** The first approach — measuring grid distortion on the
tifxyz quads — does not work, and the reason is structural rather than a matter
of tuning: when a tracer walks onto the neighbouring sheet the parameterisation
stays locally isometric, so distortion is blind to the defect by construction.
On real data it fires only along hole boundaries, which the validity mask
already gives for free. Documented in the repo so the next person does not
spend a day rediscovering it.

## Viewer

![viewer](viewer-screenshot.png)

A percentage tells you a segment is 30 % wrong. It does not tell you *where*,
or in what shape. The viewer draws the traced surfaces in 3D with escaped
regions in red, so a failure's geometry is immediately legible — in the shot
above, the whole upper scroll is red in one glance.

Meshes stream directly from the open-data bucket (it is CORS-open); only small
per-vertex flag files are served locally, so the viewer stays under 10 MB and
always shows current upstream geometry.

## Notes that may be useful regardless

- `tifxyz` is three uncompressed single-strip float32 TIFFs — a 12-line reader,
  no TIFF library.
- The volumes are uncompressed, so one z-plane inside a chunk is a contiguous
  16 KB block. A Range request pulls a full-resolution cross-section for
  ~1.3 MB instead of 162 MB.
- A missing chunk is not an error: zarr's `fill_value` is `0` here, which means
  outside the scroll. 404 is an answer.
- Meshes are registered to a specific volume (`-on-<volume-id>-<res>um.tifxyz`).
  Sampling against a different volume of the same scroll gives silent garbage.
- Scroll 1's axis wanders ~400 voxels over its height; a single global centre
  makes radius meaningless near the core.
- Winding numbers in segment names run inside-out — w001 is innermost.
  Verified: median radius rises monotonically across all 28 winding groups.

## Limitations

- Only checks whether a surface is *on the scroll*. It does not detect sheet
  switches between adjacent layers — the surface is still inside the object.
  That remains the open problem.
- Requires a masked volume, which excludes 5 scrolls that have segments.
- Absolute rates carry ~±30 % resolution uncertainty; rankings are sound.
- The distinction between "the tracer extrapolated into space" and "the mask
  excludes frayed outer layers that are really papyrus" is not resolved here.
  Both are worth knowing about; only the first is a bug.
