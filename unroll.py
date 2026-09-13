"""Split a range-named segment (w010-027) into one point set per winding.

Scroll 1 names its segments by winding range, so the adjacent-winding check has
no pairs to work with. But the range in the name is a calibration: a segment
labelled w010-027 must span exactly 18 turns. Unwrap the angle about the scroll
axis, confirm the span matches the label, and the per-vertex winding number
falls out - after which the *validated* cross-segment comparison applies, with
none of the row-slant trouble that breaks a within-row walk.
"""
import os, re, sys, glob, numpy as np
from l0 import read_tifxyz_plane
from axis import fit_axis, radius
from coverage import load_seg
from scan import CACHE

def grid(d):
    P = np.dstack([read_tifxyz_plane(f"{d}/{a}.tif") for a in "xyz"])
    v = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P, v

def range_of(name):
    m = re.search(r"w(\d{3})-(\d{3})", name)
    return (int(m.group(1)), int(m.group(2))) if m else None

def winding_field(P, valid, zc, cen):
    """Unwrapped angle in turns, consistent across rows."""
    cx = np.interp(P[..., 2], zc, cen[:, 0])
    cy = np.interp(P[..., 2], zc, cen[:, 1])
    a = np.where(valid, np.arctan2(P[..., 1] - cy, P[..., 0] - cx), np.nan)
    out = np.full(a.shape, np.nan)
    for i in range(a.shape[0]):
        m = np.isfinite(a[i])
        if m.sum() >= 3: out[i, m] = np.unwrap(a[i, m])
    turns = out / (2 * np.pi)
    # rows are unwrapped independently and can land a whole turn apart; pull each
    # onto the same branch using the column-wise median as the reference
    with np.errstate(all="ignore"):
        cols = np.isfinite(turns).any(axis=0)
        ref = np.full(turns.shape[1], np.nan)
        ref[cols] = np.nanmedian(turns[:, cols], axis=0)
    ok = np.isfinite(ref)
    for i in range(turns.shape[0]):
        m = np.isfinite(turns[i]) & ok
        if m.sum() < 3: continue
        turns[i] -= np.round(np.nanmedian(turns[i][m] - ref[m]))
    return turns

def slices(d, zc, cen, tol=0.25):
    """{winding number -> points}, or None if the span contradicts the name."""
    name = os.path.basename(d.rstrip("/"))
    rng = range_of(name)
    if not rng: return None, "no winding range in name"
    lo, hi = rng
    P, v = grid(d)
    t = winding_field(P, v, zc, cen)
    f = t[np.isfinite(t)]
    if len(f) < 500: return None, "too few vertices"
    span = np.percentile(f, 99.5) - np.percentile(f, 0.5)
    want = hi - lo + 1      # the range in the name is inclusive: w116-117 is two turns
    if want and abs(span - want) / want > tol:
        return None, f"span {span:.1f} turns vs {want} in the name"
    # inner windings carry the lower number, so orient by radius
    r = radius(P.reshape(-1, 3), zc, cen).reshape(P.shape[:2])
    inner = np.nanmean(np.where(t < np.nanpercentile(t, 20), r, np.nan))
    outer = np.nanmean(np.where(t > np.nanpercentile(t, 80), r, np.nan))
    tt = t if outer > inner else -t
    base = np.nanpercentile(tt, 0.5)
    w = lo + (tt - base)
    out = {}
    for k in range(lo, hi + 1):
        m = v & (np.floor(w) == k)      # winding k owns the turn interval [k, k+1)
        if m.sum() >= 200: out[k] = P[m]
    return out, f"span {span:.1f}/{want} turns, {len(out)} slices"

if __name__ == "__main__":
    scroll = sys.argv[1]
    dirs = [d for d in sorted(glob.glob(f"{CACHE}/{scroll}/*/"))
            if os.path.exists(f"{d}/x.tif") and range_of(os.path.basename(d.rstrip("/")))]
    allp = np.concatenate([load_seg(d) for d in dirs])
    zc, cen = fit_axis(allp)
    ok = 0
    for d in dirs:
        n = os.path.basename(d.rstrip("/"))
        sl, msg = slices(d, zc, cen)
        print(f"  {n[:40]:40s} {msg}")
        if sl: ok += 1
    print(f"\n{ok}/{len(dirs)} segments unrolled")
