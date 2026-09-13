"""Cross-section of Scroll 1's core with every traced surface drawn on top.

Answers one question: at this height, which sheets near the core has anyone
traced, and which are still bare? Pulls only the chunk neighbourhood around the
axis - the volume is uncompressed uint8, so a chunk is exactly 128^3 bytes.
"""
import os, sys, glob, urllib.request, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
from l0 import read_tifxyz_plane

BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
VOL = "PHercParis4/volumes/20260310170716-45.532um-11.0m-74keV-masked.zarr"
CH = 128
INNERMOST = "5753_-7"          # deepest segment there is - its centroid tracks the core

def load_grid(d):
    P = np.dstack([read_tifxyz_plane(f"{d}/{a}.tif") for a in "xyz"])
    valid = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P, valid

def chunk(zc, yc, xc):
    try:
        raw = urllib.request.urlopen(f"{BUCKET}/{VOL}/0/{zc}/{yc}/{xc}", timeout=120).read()
        return np.frombuffer(raw, np.uint8).reshape(CH, CH, CH)
    except Exception:
        return np.zeros((CH, CH, CH), np.uint8)      # fill_value 0 = outside the mask

def plane(z, cy, cx, half=192):
    """Assemble just enough chunks to cover a 2*half box around (cy, cx)."""
    ycs = list(range(int(cy-half)//CH, int(cy+half)//CH + 1))
    xcs = list(range(int(cx-half)//CH, int(cx+half)//CH + 1))
    jobs = [(yc, xc) for yc in ycs for xc in xcs]
    with ThreadPoolExecutor(12) as ex:
        parts = list(ex.map(lambda t: chunk(z//CH, *t), jobs))
    img = np.zeros((len(ycs)*CH, len(xcs)*CH), np.uint8)
    for (yc, xc), c in zip(jobs, parts):
        img[(ycs.index(yc))*CH:(ycs.index(yc)+1)*CH,
            (xcs.index(xc))*CH:(xcs.index(xc)+1)*CH] = c[z % CH]
    return img, ycs[0]*CH, xcs[0]*CH

def rows_near(P, valid, z, tol):
    """Grid rows whose median z sits within tol of the requested plane."""
    zs = np.where(valid, P[..., 2], np.nan)
    with np.errstate(all="ignore"):
        med = np.nanmedian(zs, axis=1)
    return np.where(np.abs(med - z) < tol)[0]

if __name__ == "__main__":
    root, z, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    grids = {os.path.basename(d.rstrip("/")): load_grid(d)
             for d in sorted(glob.glob(f"{root}/*/"))}

    key = next(k for k in grids if INNERMOST in k)
    P, v = grids[key]
    r = rows_near(P, v, z, 30)
    seed = P[r][v[r]]
    cx, cy = seed[:, 0].mean(), seed[:, 1].mean()
    print(f"z={z}  axis from {key}: ({cx:.0f}, {cy:.0f})")

    img, y0, x0 = plane(z, cy, cx)
    print(f"plane {img.shape} at y0={y0} x0={x0}; inside mask {np.count_nonzero(img)/img.size:.1%}")

    fig, ax = plt.subplots(1, 2, figsize=(17, 8.5))
    for a in ax:
        a.imshow(img, cmap="gray", origin="upper",
                 extent=[x0, x0+img.shape[1], y0+img.shape[0], y0])
        a.set_xlim(cx-170, cx+170); a.set_ylim(cy+170, cy-170)
        a.set_xlabel("x (vx)"); a.set_ylabel("y (vx)")
    ax[0].set_title(f"PHercParis4  45.5 um  z={z}   raw CT")
    n = 0
    for name, (P, v) in grids.items():
        w = "w0" in name or "w1" in name
        for i in rows_near(P, v, z, 10):
            m = v[i]
            if m.sum() < 4: continue
            ax[1].plot(P[i, m, 0], P[i, m, 1], "-", lw=1.1,
                       color="tab:orange" if w else "tab:cyan", alpha=.85)
            n += 1
    ax[1].set_title(f"traced surfaces crossing this height ({n} rows)\n"
                    "orange = w-family (w010-w129)   cyan = inner / legacy segments")
    fig.tight_layout(); fig.savefig(out, dpi=120); print("wrote", out)
