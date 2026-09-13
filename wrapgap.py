"""Sheet-switch check inside a single segment, for scrolls whose segments span
several windings.

Where segments are named by a single winding, adjacent ones can be compared
directly (sheetswitch.py). Scroll 1 names its segments by winding *range*
(w010-027), so no adjacent pair exists - but a segment covering 18 turns
contains the comparison inside itself: walk one full turn along a grid row and
you must land one papyrus thickness away, every time.

  gap collapses  -> the trace came back onto the sheet it was already on
  gap doubles    -> it skipped a wrap
  gap is steady  -> the winding structure is sound
"""
import os, sys, glob, numpy as np
from l0 import read_tifxyz_plane
from coverage import load_seg
from axis import fit_axis, radius
from scan import CACHE

def grid(d):
    P = np.dstack([read_tifxyz_plane(f"{d}/{a}.tif") for a in "xyz"])
    v = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P, v

def unwrapped_angle(P, valid, zc, cen):
    """Angle about the scroll axis, unwrapped along the wrapping direction."""
    cx = np.interp(P[..., 2], zc, cen[:, 0])
    cy = np.interp(P[..., 2], zc, cen[:, 1])
    a = np.arctan2(P[..., 1] - cy, P[..., 0] - cx)
    a = np.where(valid, a, np.nan)
    # unwrap row by row, skipping invalid cells
    out = np.full(a.shape, np.nan)
    for i in range(a.shape[0]):
        m = np.isfinite(a[i])
        if m.sum() < 3: continue
        out[i, m] = np.unwrap(a[i, m])
    return np.degrees(out)

def wrap_gaps(P, valid, ang, rad):
    """Radial distance from each vertex to the point one full turn along its row.

    Not the 3D distance: grid rows slant in z, badly so on the outer windings,
    where the landing point is millimetres away in height and the 3D distance
    measures that slant instead of the sheet spacing (2627 um vs 471 um radial
    on w098-100). Only the radial step is the inter-sheet gap.
    """
    h, w, _ = P.shape
    gap = np.full((h, w), np.nan)
    for i in range(h):
        a = ang[i]
        k = np.where(np.isfinite(a))[0]
        if len(k) < 8: continue
        aa = a[k]
        s = 1 if aa[-1] > aa[0] else -1
        x, y = (aa, k) if s > 0 else (aa[::-1], k[::-1])
        for j in k:
            t = a[j] - s * 360.0
            if t < x[0] or t > x[-1]: continue
            jf = np.interp(t, x, y)
            lo, hi = int(np.floor(jf)), int(np.ceil(jf))
            if not (valid[i, lo] and valid[i, hi]): continue
            f = jf - lo
            gap[i, j] = abs(rad[i, j] - (rad[i, lo] * (1 - f) + rad[i, hi] * f))
    return gap

def run(scroll, um, pattern="*"):
    dirs = [d for d in sorted(glob.glob(f"{CACHE}/{scroll}/{pattern}/"))
            if os.path.exists(f"{d}/x.tif")]
    allp = np.concatenate([load_seg(d) for d in dirs])
    zc, cen = fit_axis(allp)
    print(f"=== {scroll}  {len(dirs)} segments  axis over z "
          f"{zc.min():.0f}..{zc.max():.0f}  ({um} um/voxel)")
    print(f"  {'segment':34s} {'cover':>6} {'p10':>6} {'p50':>6} {'p90':>6} "
          f"{'um':>6} {'lo/med':>7} {'hi/med':>7}")
    rows = []
    for d in dirs:
        n = os.path.basename(d.rstrip("/"))
        P, v = grid(d)
        ang = unwrapped_angle(P, v, zc, cen)
        rad = radius(P.reshape(-1, 3), zc, cen).reshape(P.shape[:2])
        g = wrap_gaps(P, v, ang, rad)
        f = g[np.isfinite(g)]
        if len(f) < 200:
            print(f"  {n[:34]:34s} {'--':>6}  too few wrap pairs"); continue
        q = np.percentile(f, [10, 50, 90])
        cover = np.isfinite(g).sum() / max(v.sum(), 1)
        lo, hi = q[0] / q[1], q[2] / q[1]
        flag = ""
        if lo < 0.30: flag = "  <-- collapses somewhere"
        elif hi > 2.5: flag = "  <-- doubles somewhere"
        rows.append(dict(segment=n, cover=float(cover), p50=float(q[1]),
                         um=float(q[1] * um), lo=float(lo), hi=float(hi), flag=flag.strip(" <-")))
        print(f"  {n[:34]:34s} {cover:>6.0%} {q[0]:>6.1f} {q[1]:>6.1f} {q[2]:>6.1f} "
              f"{q[1]*um:>6.0f} {lo:>7.2f} {hi:>7.2f}{flag}")
    if rows:
        med = np.median([r["p50"] for r in rows])
        print(f"  -- median wrap gap {med:.1f} vx = {med*um:.0f} um "
              f"over {len(rows)} segments")
    return rows

if __name__ == "__main__":
    import json
    plans = {p["scroll"]: p for p in json.load(open("scan_plan.json"))}
    s = sys.argv[1]
    rows = run(s, plans[s]["base_um"], sys.argv[2] if len(sys.argv) > 2 else "*")
    json.dump(rows, open(f"wrapgap_{s}.json", "w"), indent=1)
