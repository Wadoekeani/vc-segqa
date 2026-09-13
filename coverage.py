"""Where on Scroll 1 has anyone actually traced a surface?

Deliberately model-free: no umbilicus, no winding number, no quality judgement.
Only "is there mesh here or not", which is an objective fact.
"""
import os, sys, glob, numpy as np
from l0 import read_tifxyz_plane

def load_seg(d):
    P = np.dstack([read_tifxyz_plane(f"{d}/{a}.tif") for a in "xyz"])
    valid = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P[valid]                      # (n, 3) points, order no longer matters

def load_all(root):
    segs = {}
    for d in sorted(glob.glob(f"{root}/*/")):
        name = os.path.basename(d.rstrip("/"))
        try:
            pts = load_seg(d)
            if len(pts): segs[name] = pts
        except Exception as e:
            print(f"  skip {name}: {e}")
    return segs

if __name__ == "__main__":
    segs = load_all(sys.argv[1])
    allpts = np.concatenate(list(segs.values()))
    lo, hi = allpts.min(0), allpts.max(0)
    print(f"{len(segs)} segments, {len(allpts):,} vertices")
    print(f"bbox x {lo[0]:.0f}..{hi[0]:.0f}  y {lo[1]:.0f}..{hi[1]:.0f}  z {lo[2]:.0f}..{hi[2]:.0f} vx")
    print(f"\n{'segment':34s} {'pts':>7} {'z range':>16} {'z span':>7}")
    rows = []
    for n, p in segs.items():
        rows.append((p[:, 2].min(), p[:, 2].max(), len(p), n))
    for z0, z1, np_, n in sorted(rows):
        print(f"{n:34s} {np_:>7,} {z0:>7.0f}..{z1:<7.0f} {z1-z0:>7.0f}")
