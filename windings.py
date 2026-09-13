"""Winding coverage of Scroll 1, using the official wNNN-MMM segment names."""
import os, re, sys, glob, numpy as np
from coverage import load_seg
from axis import fit_axis, radius

root = sys.argv[1]
data = {}
for d in sorted(glob.glob(f"{root}/*/")):
    n = os.path.basename(d.rstrip("/"))
    try: data[n] = load_seg(d)
    except Exception as e: print("skip", n, e)
allp = np.concatenate(list(data.values()))
zc, cen = fit_axis(allp)
print(f"{len(data)} segments, {len(allp):,} vertices")
print(f"axis fitted over z {zc.min():.0f}..{zc.max():.0f} in {len(zc)} slabs; "
      f"drift x {cen[:,0].min():.0f}..{cen[:,0].max():.0f} y {cen[:,1].min():.0f}..{cen[:,1].max():.0f} vx\n")

groups = {}
for n, p in data.items():
    m = re.search(r"w(\d{3})-(\d{3})", n)
    groups.setdefault((int(m.group(1)), int(m.group(2))) if m else None, []).append(n)

print(f"{'windings':>10} {'r p5':>6} {'p50':>6} {'p95':>6}")
prev = None
for k in sorted(x for x in groups if x):
    p = np.concatenate([data[n] for n in groups[k]])
    r = radius(p, zc, cen)
    q = np.percentile(r, [5, 50, 95])
    print(f"{k[0]:>4}-{k[1]:<4} {q[0]:>6.0f} {q[1]:>6.0f} {q[2]:>6.0f}")
    assert prev is None or q[1] > prev, "winding radius must increase monotonically"
    prev = q[1]

wp = np.concatenate([data[n] for k in groups if k for n in groups[k]])
inner = np.percentile(radius(wp, zc, cen), 1)
print(f"\nw-family inner edge (p1 radius): {inner:.0f} vx")
print(f"{'unnamed segment':40s} {'r p1':>6} {'r p50':>6} {'frac inside':>12}")
for n in sorted(groups.get(None, [])):
    r = radius(data[n], zc, cen)
    print(f"{n:40s} {np.percentile(r,1):>6.0f} {np.percentile(r,50):>6.0f} "
          f"{(r < inner).mean():>11.2%}")
