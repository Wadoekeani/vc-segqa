"""Cross-resolution agreement: one segment, meshed on two volumes, compared in
the shared reference frame.

Every fine volume ships a transform.json mapping its voxels onto the 7.91 um
volume from 2023, so that volume is the common frame - and meshes cut on it need
no transform at all, which makes them the ruler. Two independent reconstructions
of the same sheet must land on the same sheet; sheets here are ~180 um apart, so
disagreement beyond ~90 um means one of the two is on the wrong layer.

Distance is point-to-surface, not point-to-vertex: the grids are ~150 um apart,
so nearest-vertex distance has a ~80 um floor that swamps the thing being
measured (see _control in __main__).
"""
import glob, json, os, sys, numpy as np, l0

CACHE = os.path.expanduser(os.environ.get("SEGQA_CACHE", "~/.cache/segqa"))
TDIR = f"{CACHE}/_transforms"
UM = 7.91                                   # fixed-frame voxel size
FIXED_VOL = "20230205180739"
VOLS = {"45.5": "20260310170716-45.532um-11.0m-74keV-masked",
        "2.4": "20260411134726-2.400um-0.2m-78keV-masked",
        "1.129": "20260608103018-1.129um-0.2m-78keV-masked"}


def find_mesh(seg, res):
    """Cache dir for one segment at one resolution, or None.

    Coarse meshes live where the rest of segqa puts them; the finer ones are
    under _fine/ with the mesh folder's own name, which is not always prefixed
    with the segment id - so match on substring rather than rebuild the name.
    """
    if res == "45.5":
        for d in (f"{CACHE}/PHercParis4/{seg}", f"{CACHE}/meshes/{seg}"):
            if os.path.exists(f"{d}/x.tif"): return d
        return None
    hits = [d for d in sorted(glob.glob(f"{CACHE}/_fine/*{res}um.tifxyz"))
            if seg in os.path.basename(d) and os.path.exists(f"{d}/x.tif")]
    return hits[0] if hits else None


def to_fixed(vol):
    """voxel coords of `vol` -> fixed-frame coords. Direction is M @ moving,
    which the landmark pairs in the file confirm (the other way is off by 40k)."""
    M = np.array(json.load(open(f"{TDIR}/{vol}.json"))["transformation_matrix"])
    def f(P):
        with np.errstate(all="ignore"):      # Accelerate BLAS raises spurious flags
            return np.asarray(P, np.float64) @ M[:, :3].T + M[:, 3]
    return f


def _normal(P, valid, i, j):
    """Sheet normal at grid node (i,j), oriented by the grid's own parameter
    directions so the sign of a distance is comparable between neighbours - an
    SVD normal flips arbitrarily and makes "is this an offset?" unanswerable."""
    H, W = valid.shape
    if not (0 < i < H-1 and 0 < j < W-1): return None
    if not (valid[i, j-1] and valid[i, j+1] and valid[i-1, j] and valid[i+1, j]): return None
    du = P[i, j+1].astype(np.float64) - P[i, j-1]
    dv = P[i+1, j].astype(np.float64) - P[i-1, j]
    n = np.cross(du, dv); L = np.linalg.norm(n)
    return None if L < 1e-9 else n / L


def _buckets(P, valid, cell, lo):
    b = {}
    ij = np.argwhere(valid)
    for i, j in ij:
        k = tuple(np.floor((P[i, j] - lo) / cell).astype(np.int64))
        b.setdefault(k, []).append((i, j))
    return b


def surface_dist(Q, P, valid, cell=40.0):
    """Distance from each point of Q to the sheet P (grid + valid mask).

    Nearest vertex first, then the distance to a plane fitted through that
    vertex's 3x3 grid neighbourhood - otherwise the grid spacing is the floor.
    """
    lo = P[valid].min(0).astype(np.float64)
    b = _buckets(P, valid, cell, lo)
    off = [(a, c, e) for a in (-1, 0, 1) for c in (-1, 0, 1) for e in (-1, 0, 1)]
    H, W = valid.shape
    out = np.full(len(Q), np.nan)
    for n, q in enumerate(Q):
        k = tuple(np.floor((q - lo) / cell).astype(np.int64))
        cand = [p for a, c, e in off for p in b.get((k[0]+a, k[1]+c, k[2]+e), ())]
        if not cand: continue
        cand = np.array(cand)
        pts = P[cand[:, 0], cand[:, 1]].astype(np.float64)
        i, j = cand[np.linalg.norm(pts - q, axis=1).argmin()]
        i0, i1, j0, j1 = max(0, i-1), min(H, i+2), max(0, j-1), min(W, j+2)
        patch = P[i0:i1, j0:j1][valid[i0:i1, j0:j1]].astype(np.float64)
        if len(patch) < 4:                       # edge of the sheet: vertex distance
            out[n] = np.linalg.norm(P[i, j].astype(np.float64) - q); continue
        c = patch.mean(0)
        nrm = _normal(P, valid, i, j)
        if nrm is None:
            nrm = np.linalg.svd(patch - c)[2][2]  # edge: sign is then arbitrary
        out[n] = (q - c) @ nrm
    return out


def load_fixed(path, res):
    """A mesh in fixed-frame coords. res '7.91' is already in that frame."""
    P, valid = l0.load(path)
    if res == "7.91":
        return np.asarray(P, np.float64), valid
    return to_fixed(VOLS[res])(P.reshape(-1, 3)).reshape(P.shape), valid


def compare(patha, ra, pathb, rb, n=3000):
    Pa, va = load_fixed(patha, ra)
    Pb, vb = load_fixed(pathb, rb)
    Q = Pa[va]; Q = Q[::max(1, len(Q)//n)]
    d = surface_dist(Q, Pb, vb)
    return Q[~np.isnan(d)], d[~np.isnan(d)]


if __name__ == "__main__":
    seg = sys.argv[1] if len(sys.argv) > 1 else "20230702185753"
    p45, p79, p24 = (find_mesh(seg, r) for r in ("45.5", "7.91", "2.4"))
    assert p79, f"no 7.91 um mesh cached for {seg} (run xres_fetch.py 7.91um)"
    P, v = load_fixed(p79, "7.91")
    # control: midpoints between grid nodes lie on the sheet by construction, so
    # whatever this reports is the measurement's own noise, not disagreement
    m = v[:, :-1] & v[:, 1:]
    mid = (P[:, :-1][m] + P[:, 1:][m]) / 2
    c = surface_dist(mid[::max(1, len(mid)//2000)], P, v)
    c = np.abs(c[~np.isnan(c)]) * UM
    print(f"control (on-sheet midpoints vs its own mesh): p50={np.median(c):.1f} "
          f"p90={np.percentile(c,90):.1f} max={c.max():.1f} um")
    assert np.median(c) < 10, f"measurement floor too high ({np.median(c):.1f} um)"

    print(f"\n{seg}   (um, and voxels of the coarser volume in the pair)")
    for (a, pa), (b, pb) in [(("45.5", p45), ("7.91", p79)), (("2.4", p24), ("7.91", p79)),
                             (("45.5", p45), ("2.4", p24))]:
        if not (pa and pb):
            print(f"  {a:>5} vs {b:<5} skipped: mesh not cached"); continue
        Q, sd = compare(pa, a, pb, b)
        d = np.abs(sd) * UM
        VX = {"45.5": 45.532, "2.4": 2.4, "7.91": 7.91}
        coarse = a if VX[a] >= VX[b] else b        # the coarser of the two sets the scale
        print(f"  {a:>5} vs {b:<5} n={len(d):<5} p50={np.median(d):6.1f} "
              f"p90={np.percentile(d,90):6.1f} max={d.max():7.1f} um"
              f"   p50 = {np.median(d)/VX[coarse]:.2f} x {coarse} um voxel")
        # offset or wrong shape? block the points and look at signed distance:
        # a big mean with a small spread is a local shift, the reverse is shape
        blk = np.floor(Q / 1000).astype(np.int64)      # ~8 mm cells
        keys, inv = np.unique(blk, axis=0, return_inverse=True)
        mu = np.array([sd[inv == k].mean() for k in range(len(keys))]) * UM
        sg = np.array([sd[inv == k].std() for k in range(len(keys))]) * UM
        cnt = np.bincount(inv)
        big = cnt >= 25
        if big.sum():
            print(f"         {big.sum()} blocks (>=25 pts): |mean| p50={np.median(abs(mu[big])):.1f} "
                  f"p90={np.percentile(abs(mu[big]),90):.1f} | within-block sd p50={np.median(sg[big]):.1f} um")
