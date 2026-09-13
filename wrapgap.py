"""Inter-wrap spacing: umbilicus-free sheet-switch detector.

Walk one full turn along the grid row. The point you land on must sit on the
NEXT sheet, one papyrus thickness away. If the tracer switched sheets, that gap
collapses towards zero (landed on the same sheet) or roughly doubles (skipped one).
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from l0 import load

def wrap_gaps(P, valid):
    h, w, _ = P.shape
    xy = P[..., :2]
    # per-row centre only to define an angle; the gap itself is measured in 3D
    C = np.full((h, 2), np.nan)
    for i in range(h):
        p = xy[i][valid[i]]
        if len(p) >= 8:
            A = np.c_[2*p, np.ones(len(p))]
            C[i] = np.linalg.lstsq(A, (p**2).sum(1), rcond=None)[0][:2]
    ang = np.degrees(np.unwrap(np.arctan2(xy[..., 1]-C[:, None, 1],
                                          xy[..., 0]-C[:, None, 0]), axis=1))
    gap = np.full((h, w), np.nan)
    for i in range(h):
        if not np.isfinite(C[i]).all(): continue
        a, v = ang[i], valid[i]
        s = -1 if a[v].argmin() > a[v].argmax() else 1   # winding direction
        for j in range(w):
            if not v[j]: continue
            target = a[j] - s*360.0
            k = np.where(v)[0]
            aa = a[k]
            if target > aa.max() or target < aa.min(): continue
            j2 = np.interp(target, aa[::s] if s < 0 else aa, k[::s] if s < 0 else k)
            lo, hi = int(np.floor(j2)), int(np.ceil(j2))
            if not (v[lo] and v[hi]): continue
            t = j2 - lo
            Q = P[i, lo]*(1-t) + P[i, hi]*t
            gap[i, j] = np.linalg.norm(P[i, j] - Q)
    return gap

def main(path, out, um):
    P, valid = load(path)
    gap = wrap_gaps(P, valid)
    g = gap[np.isfinite(gap)]
    if len(g) < 50: print("too few wrap pairs"); return
    med = np.median(g)
    print(f"wrap pairs: {len(g)}  ({np.isfinite(gap).mean():.1%} of grid)")
    print(f"inter-wrap gap: p5={np.percentile(g,5):.1f} p50={med:.1f} "
          f"p95={np.percentile(g,95):.1f} vx"
          f"   -> p50 = {med*um:.0f} um")
    rel = gap/med
    for name, sel in [("collapsed (<0.4x)", rel < 0.4), ("doubled (>1.7x)", rel > 1.7)]:
        n = np.nansum(sel)
        print(f"  {name:20s} {int(n):5d} quads ({n/len(g):.2%})")

    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    ax[0].hist(g, bins=80); ax[0].axvline(med, color="r")
    ax[0].set_xlabel("inter-wrap gap (vx)"); ax[0].set_title("gap distribution")
    im = ax[1].imshow(rel, cmap="coolwarm", vmin=0, vmax=2,
                      interpolation="nearest", aspect="auto")
    ax[1].set_title("gap / median"); fig.colorbar(im, ax=ax[1], shrink=.7)
    fig.tight_layout(); fig.savefig(out, dpi=110); print("wrote", out)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 45.532)
