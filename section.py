"""Render any z cross-section of Scroll 1, with traced surfaces overlaid.

The volume is uncompressed uint8 in C order, so one z-plane inside a 128^3
chunk is a contiguous 16 KB block: fetch it with a Range request instead of
pulling the whole 2 MB chunk. A full native-resolution slice costs ~5 MB.
"""
import os, sys, glob, urllib.request, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
from l0 import read_tifxyz_plane

BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
VOL = "PHercParis4/volumes/20260310170716-45.532um-11.0m-74keV-masked.zarr"
CH = 128
SHAPE = {0: 2264, 1: 1132, 2: 566, 3: 283, 4: 142, 5: 71}   # y=x extent per level

def fetch_plane_chunk(level, zc, zl, yc, xc):
    off = zl * CH * CH
    req = urllib.request.Request(f"{BUCKET}/{VOL}/{level}/{zc}/{yc}/{xc}",
                                 headers={"Range": f"bytes={off}-{off+CH*CH-1}"})
    try:
        return np.frombuffer(urllib.request.urlopen(req, timeout=60).read(),
                             np.uint8).reshape(CH, CH)
    except Exception:
        return np.zeros((CH, CH), np.uint8)     # fill_value 0 = outside the mask

def plane(level, z):
    """Full z-plane at the given pyramid level, as (img, n_bytes)."""
    zl_full = z >> level
    zc, zl = divmod(zl_full, CH)
    n = (SHAPE[level] + CH - 1) // CH
    jobs = [(yc, xc) for yc in range(n) for xc in range(n)]
    with ThreadPoolExecutor(16) as ex:
        parts = list(ex.map(lambda t: fetch_plane_chunk(level, zc, zl, *t), jobs))
    img = np.zeros((n*CH, n*CH), np.uint8)
    for (yc, xc), p in zip(jobs, parts):
        img[yc*CH:(yc+1)*CH, xc*CH:(xc+1)*CH] = p
    return img[:SHAPE[level], :SHAPE[level]], len(jobs) * CH * CH

def load_grid(d):
    P = np.dstack([read_tifxyz_plane(f"{d}/{a}.tif") for a in "xyz"])
    v = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P, v

def rows_near(P, v, z, tol):
    zs = np.where(v, P[..., 2], np.nan)
    rows = v.any(axis=1)
    med = np.full(len(zs), np.inf)
    med[rows] = np.nanmedian(zs[rows], axis=1)
    return np.where(np.abs(med - z) < tol)[0]

def draw(ax, grids, z, level, tol=10):
    n = 0
    for name, (P, v) in grids.items():
        w = "-w0" in name or "-w1" in name
        for i in rows_near(P, v, z, tol):
            m = v[i]
            if m.sum() < 4: continue
            ax.plot(P[i, m, 0] / (1 << level), P[i, m, 1] / (1 << level), "-",
                    lw=.8, color="tab:orange" if w else "tab:cyan", alpha=.9)
            n += 1
    return n

if __name__ == "__main__":
    root, level, out = sys.argv[1], int(sys.argv[2]), sys.argv[-1]
    zs = [int(a) for a in sys.argv[3:-1]]
    grids = {os.path.basename(d.rstrip("/")): load_grid(d)
             for d in sorted(glob.glob(f"{root}/*/"))}
    fig, ax = plt.subplots(2, len(zs), figsize=(6.2*len(zs), 12.4), squeeze=False)
    for j, z in enumerate(zs):
        img, nb = plane(level, z)
        print(f"z={z:5d} level={level}  plane {img.shape}  {nb/1e6:.1f} MB  "
              f"inside mask {np.count_nonzero(img)/img.size:.1%}")
        for i in (0, 1):
            ax[i][j].imshow(img, cmap="gray", origin="upper")
            ax[i][j].set_xticks([]); ax[i][j].set_yticks([])
        ax[0][j].set_title(f"z = {z}   raw CT")
        n = draw(ax[1][j], grids, z, level)
        ax[1][j].set_title(f"z = {z}   {n} traced rows")
    fig.tight_layout(); fig.savefig(out, dpi=110); print("wrote", out)
