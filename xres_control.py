"""Control arm for xres_scan.py: the fine volumes compared against each other.

xres_scan reports that a segment's 45.5 um mesh and its 7.91 um mesh disagree by
about one 45.5 um voxel. That number only means something if the comparison
itself is tighter than that - the transform could be wrong, the point-to-surface
distance could be wrong, the two meshes might not even be the same sheet.

So run the same comparison where the answer should be small: 2.4 vs 7.91, and
1.129 vs 2.4. Both are the same kind of pair (one segment, two volumes, one
official transform) and neither involves the coarse volume.

The 1.129 um transform targets the 2.4 um volume rather than the 7.91 um one, so
that pair is compared in the 2.4 volume's own frame.
"""
import json, os, numpy as np, l0, xres
C = xres.CACHE; F = f"{C}/_fine"

def load_in(path, vol=None, n=None):
    P, v = l0.load(path)
    if vol is None: return np.asarray(P, np.float64), v
    return xres.to_fixed(vol)(P.reshape(-1, 3)).reshape(P.shape), v

def cmp(Pa, va, Pb, vb, n=2500):
    Q = Pa[va]; Q = Q[::max(1, len(Q)//n)]
    d = xres.surface_dist(Q, Pb, vb)
    return np.abs(d[~np.isnan(d)])

# One per group of the xres_scan breakdown, so the control is not accidentally
# sampled only from the group that scores best at 45.5 um.
SEGS = [("20230702185753", "2023"), ("20230929220926", "2023"),
        ("20231005123336", "2023"), ("20231031143852", "2023"),
        ("20231210121321", "2023"),
        ("20260602204401-5753_-7", "5753"), ("20260603024952-5753_-2", "5753"),
        ("20260623141924-w010-027", "w-named, inner"),
        ("20260623171929-w128-129", "w-named, outer")]

for seg, label in SEGS:
    p79, p24 = xres.find_mesh(seg, "7.91"), xres.find_mesh(seg, "2.4")
    if not (p79 and p24):
        print(f"{seg}: missing a mesh (run xres_fetch.py 2.4um)"); continue
    A, va = load_in(p24, xres.VOLS["2.4"])
    B, vb = load_in(p79)
    d = cmp(A, va, B, vb) * xres.UM
    print(f"{seg:<24} {label:<15} 2.4 vs 7.91: n={len(d):<5} p50={np.median(d):6.2f} "
          f"p90={np.percentile(d,90):6.2f} p99={np.percentile(d,99):7.2f} um "
          f"= {np.median(d)/7.91:.2f} x 7.91um voxel", flush=True)

# third point on the ladder: 1.129 vs 2.4, done in the 2.4 volume's own frame
seg = "20230702185753"
p11, p24 = xres.find_mesh(seg, "1.129"), xres.find_mesh(seg, "2.4")
if p11 and p24:
    A, va = load_in(p11, xres.VOLS["1.129"])      # its transform targets the 2.4 volume
    B, vb = load_in(p24)
    d = cmp(A, va, B, vb) * 2.4
    print(f"{seg}  1.129 vs 2.4: n={len(d):<5} p50={np.median(d):6.2f} p90={np.percentile(d,90):6.2f} "
          f"p99={np.percentile(d,99):7.2f} um   = {np.median(d)/2.4:.2f} x 2.4um voxel")
