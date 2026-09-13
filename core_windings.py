"""Which absolute windings do Scroll 1's core segments actually cover?

The w-family carries official winding numbers; the 5753 family and the 2023
Grand-Prize segments do not - their names are offsets from a seed. Unrolling the
w-family gives a labelled point cloud to measure against: for every vertex of an
unlabelled segment, find the nearest labelled vertex and read off its winding.

Distance says how to read the match: about zero means the same sheet, about one
sheet spacing means the neighbouring winding - which is how a segment lying
*inside* w010 reveals itself.
"""
import os, re, sys, glob, json, numpy as np
from collections import Counter
from unroll import slices, range_of
from neighbour import nn_dist
from coverage import load_seg
from axis import fit_axis, radius
from scan import CACHE

UM = 45.532
SCROLL = "PHercParis4"

def build_reference():
    dirs = [d for d in sorted(glob.glob(f"{CACHE}/{SCROLL}/*/"))
            if os.path.exists(f"{d}/x.tif") and range_of(os.path.basename(d.rstrip("/")))]
    allp = np.concatenate([load_seg(d) for d in dirs])
    zc, cen = fit_axis(allp)
    pts, lab = [], []
    for d in dirs:
        sl, _ = slices(d, zc, cen)
        if not sl: continue
        for k, p in sl.items():
            pts.append(p); lab.append(np.full(len(p), k))
    return np.concatenate(pts), np.concatenate(lab), zc, cen

if __name__ == "__main__":
    ref, lab, zc, cen = build_reference()
    print(f"reference: {len(ref):,} labelled vertices, windings "
          f"{lab.min()}–{lab.max()}")
    inner_edge = np.percentile(radius(ref[lab == lab.min()], zc, cen), 1)
    print(f"innermost labelled winding w{lab.min():03d}: inner edge {inner_edge:.0f} vx "
          f"= {inner_edge*UM/1000:.1f} mm from the axis\n")

    rng = np.random.default_rng(0)
    others = [d for d in sorted(glob.glob(f"{CACHE}/{SCROLL}/*/"))
              if os.path.exists(f"{d}/x.tif")
              and not range_of(os.path.basename(d.rstrip("/")))]
    print(f"{'segment':40s} {'r p50':>7} {'match':>7} {'gap vx':>7} {'reading':>22}")
    rows = []
    for d in others:
        n = os.path.basename(d.rstrip("/"))
        p = load_seg(d)
        sub = p[rng.choice(len(p), min(2500, len(p)), replace=False)]
        dist, who = nn_dist(sub, ref, cell=max(8.0, 400.0 / UM), want_index=True)
        ok = np.isfinite(dist) & (who >= 0)
        if ok.sum() < 50:
            print(f"{n[:40]:40s} {'--':>7}  no labelled neighbour in range"); continue
        w = Counter(lab[who[ok]]).most_common(1)[0][0]
        gap = float(np.median(dist[ok]))
        r50 = float(np.median(radius(p, zc, cen)))
        rin = float(np.median(radius(ref[lab == w], zc, cen)))
        side = "inside" if r50 < rin else "outside"
        note = (f"= w{w:03d}" if gap < 3 else
                f"{side} w{w:03d} by {gap*UM:.0f} um")
        print(f"{n[:40]:40s} {r50:>7.0f} {w:>7} {gap:>7.1f} {note:>22}")
        rows.append(dict(segment=n, match=int(w), gap_vx=gap, gap_um=gap*UM,
                         r50=r50, side=side))
    json.dump(rows, open("core_windings.json", "w"), indent=1)
