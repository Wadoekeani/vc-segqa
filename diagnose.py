"""Diagnose flagged quads: real defect, or artifact of my metric near holes?"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import deque
from l0 import load, l0

def blobs(mask):
    lab = np.zeros(mask.shape, int); n = 0; seen = np.zeros(mask.shape, bool)
    for y, x in zip(*np.nonzero(mask)):
        if seen[y, x]: continue
        n += 1; q = deque([(y, x)]); seen[y, x] = True
        while q:
            cy, cx = q.popleft(); lab[cy, cx] = n
            for dy, dx in ((1,0),(-1,0),(0,1),(0,-1)):
                ny, nx = cy+dy, cx+dx
                if 0 <= ny < mask.shape[0] and 0 <= nx < mask.shape[1] \
                   and mask[ny, nx] and not seen[ny, nx]:
                    seen[ny, nx] = True; q.append((ny, nx))
    return lab, n

def dist_to_false(valid):
    """Chebyshev distance from every cell to the nearest invalid cell (BFS)."""
    d = np.full(valid.shape, 1 << 20, int)
    q = deque()
    for y, x in zip(*np.nonzero(~valid)):
        d[y, x] = 0; q.append((y, x))
    # treat the grid border as invalid too - holes and edges are the same story
    for y in range(valid.shape[0]):
        for x in (0, valid.shape[1]-1):
            if d[y, x] > 0: d[y, x] = 0; q.append((y, x))
    for x in range(valid.shape[1]):
        for y in (0, valid.shape[0]-1):
            if d[y, x] > 0: d[y, x] = 0; q.append((y, x))
    while q:
        cy, cx = q.popleft()
        for dy in (-1,0,1):
            for dx in (-1,0,1):
                ny, nx = cy+dy, cx+dx
                if 0 <= ny < d.shape[0] and 0 <= nx < d.shape[1] and d[ny, nx] > d[cy, cx]+1:
                    d[ny, nx] = d[cy, cx]+1; q.append((ny, nx))
    return d

def main(path, out):
    P, valid = load(path)
    m = l0(P, valid)
    s = np.nan_to_num(m["stretch"], nan=0.0)
    bad = s > 5
    lab, n = blobs(bad)
    dist = dist_to_false(valid)[:-1, :-1]        # align to quad grid
    med = m["median_edge_vx"]

    print(f"grid {m['grid']}  valid {m['valid_frac']:.1%}  median edge {med:.2f} vx")
    print(f"flagged {bad.sum()} quads in {n} blobs\n")
    print(f"{'blob':>4} {'size':>4} {'row':>5} {'col':>4} {'maxStretch':>10} "
          f"{'jump_vx':>8} {'d_to_hole':>9}")
    interior = 0
    for i in range(1, n+1):
        sel = lab == i
        ys, xs = np.nonzero(sel)
        dmin = dist[sel].min()
        if dmin >= 2: interior += 1
        print(f"{i:>4} {sel.sum():>4} {int(ys.mean()):>5} {int(xs.mean()):>4} "
              f"{s[sel].max():>10.1f} {s[sel].max()*med:>8.0f} {dmin:>9}")
    print(f"\ninterior blobs (>=2 cells from any hole/border): {interior}/{n}")

    fig, ax = plt.subplots(1, 3, figsize=(13, 7))
    ax[0].imshow(valid, cmap="gray", interpolation="nearest")
    ax[0].set_title(f"valid mask ({m['valid_frac']:.0%})")
    im = ax[1].imshow(np.where(s > 0, np.log10(s), np.nan), cmap="inferno",
                      interpolation="nearest", vmin=0, vmax=1.6)
    ax[1].set_title("log10 stretch"); fig.colorbar(im, ax=ax[1], shrink=.5)
    ax[2].imshow(valid[:-1, :-1], cmap="gray", interpolation="nearest", alpha=.5)
    ax[2].imshow(np.where(bad, lab, np.nan), cmap="tab10", interpolation="nearest")
    ax[2].set_title(f"{n} flagged blobs")
    for a in ax: a.set_xlabel("grid u"); a.set_ylabel("grid v")
    fig.tight_layout(); fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "inspect.png")
