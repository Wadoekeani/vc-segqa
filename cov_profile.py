"""Coverage along the scroll axis: how deep does tracing reach at each height?"""
import os, sys, glob, numpy as np
from coverage import load_seg
from axis import fit_axis, radius

root = sys.argv[1]
segs = {os.path.basename(d.rstrip("/")): load_seg(d)
        for d in sorted(glob.glob(f"{root}/*/"))}
allp = np.concatenate(list(segs.values()))
zc, cen = fit_axis(allp)
r = radius(allp, zc, cen)

step = 200
edges = np.arange(allp[:, 2].min(), allp[:, 2].max() + step, step)
print(f"{'z range':>14} {'segs':>5} {'verts':>9} {'r min':>6} {'r p50':>6} {'r max':>6}")
for a, b in zip(edges[:-1], edges[1:]):
    m = (allp[:, 2] >= a) & (allp[:, 2] < b)
    if m.sum() < 100:
        print(f"{a:>6.0f}..{b:<6.0f} {'-':>5} {m.sum():>9,}")
        continue
    ns = sum(1 for p in segs.values() if ((p[:, 2] >= a) & (p[:, 2] < b)).sum() > 50)
    rr = r[m]
    print(f"{a:>6.0f}..{b:<6.0f} {ns:>5} {m.sum():>9,} "
          f"{np.percentile(rr,1):>6.0f} {np.percentile(rr,50):>6.0f} {np.percentile(rr,99):>6.0f}")
