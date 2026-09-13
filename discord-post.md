Automating segmentation QA runs into a wall: validating a metric needs segments
already labelled good or bad, and there aren't any — so every metric ends up
unfalsifiable. There's a way around it. The published volumes are **masked** — voxels outside
the scroll body are exactly 0. So "does this vertex land on the object at all?"
is answerable with no labels, and a surface that leaves the scroll is wrong
by inspection.

Ran it over every published segment with a masked volume: **299 segments, 11
scrolls, 147M vertices**. Sharply bimodal — four scrolls sit at 0.00–0.25%,
PHerc1447 at 33% and PHerc1667 at 19% (17 of its 20 segments over 5%). That
four can score ~0.00% is the argument this measures something real rather than
boundary noise.

A second check, for what the first is blind to: a tracer that slips onto the
neighbouring sheet stays inside the scroll, so the mask says nothing. But
segment names carry winding numbers, and two segments one winding apart must
sit one papyrus thickness apart. That expectation is free too.

**It flags PHerc0139 w46 as having jumped onto w45's sheet:**
```
w44/45  1.28  normal
w45/46  0.25  <- coincident
w46/47  2.09  <- double gap
w47/48  1.11  normal
```
Self-consistent, and it names the culprit: had w45 moved, w44/45 would read
double, and it doesn't. Rendered, the two surfaces interpenetrate across
their whole extent. Would welcome a sanity check from anyone who knows that
scroll — real defect, or an artefact of how I'm measuring?

The two checks turn out to be independent: PHerc1667 is worst by mask escape
yet perfectly clean on windings; PHerc0139 is among the cleanest by mask and is
the only one with a sheet switch. Neither suffices alone.

Scroll 1's winding topology is also clean across 230 adjacent pairs, which
rules out mis-tracing as the explanation for the missing title there.

Code, viewer, and the two dead ends I hit (grid distortion is structurally
blind to sheet switches): https://github.com/Wadoekeani/vc-segqa
