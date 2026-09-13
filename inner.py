"""Does ANY segment reach inside the innermost traced winding (w010)?"""
import os, re, sys, glob, numpy as np
from coverage import load_seg

root = sys.argv[1]
data = {}
for d in sorted(glob.glob(f"{root}/*/")):
    n = os.path.basename(d.rstrip("/"))
    try: data[n] = load_seg(d)
    except Exception as e: print("skip", n, e)

allp = np.concatenate(list(data.values()))
cx, cy = allp[:, 0].mean(), allp[:, 1].mean()
def rad(p): return np.hypot(p[:, 0]-cx, p[:, 1]-cy)

w010 = np.concatenate([p for n, p in data.items() if "w010-027" in n])
inner_ref = np.percentile(rad(w010), 1)
print(f"centre ({cx:.0f},{cy:.0f});  w010-027 inner edge (p1 radius) = {inner_ref:.0f} vx\n")
print(f"{'segment':36s} {'pts':>8} {'r p1':>6} {'r p50':>6} {'pts inside w010':>16}")
hits = []
for n, p in sorted(data.items()):
    if re.search(r"w\d{3}-\d{3}", n): continue          # the w-family is already mapped
    r = rad(p)
    inside = int((r < inner_ref).sum())
    print(f"{n:36s} {len(p):>8,} {np.percentile(r,1):>6.0f} {np.percentile(r,50):>6.0f} "
          f"{inside:>10,} ({inside/len(p):>5.1%})")
    if inside: hits.append((n, inside))
print(f"\nsegments with ANY vertex inside w010's inner edge: {len(hits)}")
for n, c in sorted(hits, key=lambda t: -t[1])[:10]: print(f"  {n}: {c:,}")
