Automating segmentation QA runs into a wall: validating a metric needs segments
already labelled good or bad, and there aren't any — so every metric ends up
unfalsifiable. The way around it: only ask questions the published data already
answers.

**1. The volumes are masked** — outside the scroll body is exactly 0, so "is
this vertex on the object?" needs no labels. Over **299 segments, 11 scrolls,
147M vertices**: sharply bimodal, four scrolls at 0.00–0.25% against PHerc1447
at 33% and PHerc1667 at 19%.

**2. Names carry winding numbers** — two segments one winding apart must sit one
papyrus thickness apart. It flags PHerc0139 w46 as having jumped onto w45's
sheet:
```
w44/45  1.28  normal
w45/46  0.25  <- coincident
w46/47  2.09  <- double gap
w47/48  1.11  normal
```
Self-consistent, and it names the culprit: had w45 moved, w44/45 would read
double. Rendered, the two surfaces interpenetrate. Would welcome a sanity
check from anyone who knows that scroll.

**3. Each fine volume ships a transform.json** — so one segment meshed on two
volumes must describe the same sheet in the shared frame. Across 55 Scroll 1
segments, 45.5 µm vs 7.91 µm disagree by 46.6 µm median: one coarse voxel.
Controls: 2.4 vs 7.91 gives 3.7 µm, 1.129 vs 2.4 gives 0.34 µm. **Nine `5753_*`
segments sit at 82–90 µm while the other 46 sit at 33–64** — no overlap, a whole
batch flagged without anyone labelling a vertex. Two of the nine come back
among the tightest at 2.4 µm, so it's their 45.5 µm mesh, not the segments.

The three are independent: PHerc1667 is worst by mask yet perfectly clean on
windings, PHerc0139 is among the cleanest by mask and is the only sheet switch,
and check 3 catches a batch neither of the others sees. Scroll 1's topology is
clean across 230 adjacent pairs, which rules out mis-tracing as the explanation
for the missing title there.

Code, viewer, and the dead ends I kept: https://github.com/Wadoekeani/vc-segqa
