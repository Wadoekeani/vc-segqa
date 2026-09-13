"""Is there untraced papyrus inside the innermost traced winding of Scroll 1?

w001-w009 carry no segment. That is only interesting if there is material
there: the centre of a scroll can equally be a void where the umbilicus was.
Sample the masked CT inside the innermost traced surface and see.
"""
import sys, glob, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from unroll import slices, range_of
from coverage import load_seg
from axis import fit_axis, radius
from scan import CACHE
from section import plane

UM, SCROLL = 45.532, "PHercParis4"

def innermost(zc, cen):
    """Points of the lowest labelled winding, per z."""
    best = None
    for d in sorted(glob.glob(f"{CACHE}/{SCROLL}/*/")):
        n = os.path.basename(d.rstrip("/"))
        if not range_of(n) or not os.path.exists(f"{d}/x.tif"): continue
        if not n.startswith("20260701"): continue          # one batch, denser
        sl, _ = slices(d, zc, cen)
        if not sl: continue
        k = min(sl)
        if best is None or k < best[0]: best = (k, sl[k])
    return best

if __name__ == "__main__":
    dirs = [d for d in sorted(glob.glob(f"{CACHE}/{SCROLL}/*/"))
            if os.path.exists(f"{d}/x.tif") and range_of(os.path.basename(d.rstrip("/")))]
    allp = np.concatenate([load_seg(d) for d in dirs])
    zc, cen = fit_axis(allp)
    k, pts = innermost(zc, cen)
    print(f"innermost traced winding: w{k:03d}  ({len(pts):,} vertices)")
    # w010 is a thin strip, so build a radius-vs-height profile once and
    # interpolate rather than slicing it thinly per z
    r_all = radius(pts, zc, cen)
    zb = np.arange(pts[:, 2].min(), pts[:, 2].max() + 100, 100)
    idx = np.clip(np.searchsorted(zb, pts[:, 2]) - 1, 0, len(zb) - 2)
    prof_z, prof_r = [], []
    for b in range(len(zb) - 1):
        m = idx == b
        if m.sum() >= 8:
            prof_z.append((zb[b] + zb[b + 1]) / 2)
            prof_r.append(np.percentile(r_all[m], 25))
    prof_z, prof_r = np.array(prof_z), np.array(prof_r)
    print(f"  radius profile over z {prof_z.min():.0f}..{prof_z.max():.0f} "
          f"({len(prof_z)} bins), r25 {prof_r.min():.0f}..{prof_r.max():.0f} vx")

    zs = [int(z) for z in (1200, 2000, 2800)]
    fig, ax = plt.subplots(1, len(zs), figsize=(6 * len(zs), 6.4))
    for j, z in enumerate(zs):
        ring = pts[np.abs(pts[:, 2] - z) < 150]
        cx = float(np.interp(z, zc, cen[:, 0])); cy = float(np.interp(z, zc, cen[:, 1]))
        rin = float(np.interp(z, prof_z, prof_r))    # conservative inner radius
        img, _ = plane(1, z)                         # 91 um level
        L = 1
        yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
        rr = np.hypot(xx - cx / (1 << L), yy - cy / (1 << L)) * (1 << L)
        inside = rr < rin
        filled = (img > 0) & inside
        frac = filled.sum() / max(inside.sum(), 1)
        area_mm2 = filled.sum() * (UM * (1 << L) / 1000) ** 2
        print(f"  z={z}: w{k:03d} inner radius {rin:.0f} vx = {rin*UM/1000:.1f} mm; "
              f"inside it {frac:.0%} of the disc is masked-in material "
              f"({area_mm2:.1f} mm2)")
        a = ax[j]
        a.imshow(img, cmap="gray", origin="upper")
        a.plot(ring[:, 0] / (1 << L), ring[:, 1] / (1 << L), ".", ms=1.2,
               color="tab:cyan", label=f"w{k:03d} (innermost traced)")
        th = np.linspace(0, 2 * np.pi, 200)
        a.plot((cx + rin * np.cos(th)) / (1 << L), (cy + rin * np.sin(th)) / (1 << L),
               "-", lw=1, color="tab:red", label="inner radius")
        a.set_xlim((cx - 220) / (1 << L), (cx + 220) / (1 << L))
        a.set_ylim((cy + 220) / (1 << L), (cy - 220) / (1 << L))
        a.set_xticks([]); a.set_yticks([])
        a.set_title(f"z = {z}   {frac:.0%} material inside w{k:03d}")
        if j == 0: a.legend(loc="lower left", fontsize=8)
    fig.tight_layout(); fig.savefig("core-gap.png", dpi=120)
    print("wrote core-gap.png")
