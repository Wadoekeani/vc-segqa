"""Every Scroll 1 segment that exists on two volumes: do the two agree?

Writes xres_results.json. The per-segment centroid is recorded because it
separates the two candidate causes - a transform that is globally too stiff
makes the error a function of where you are in the scroll, a per-segment
meshing problem does not.

Disagreement is reported in um and in voxels of the coarse volume, and
deliberately NOT as a "fraction beyond half a sheet": every sheet-spacing number
available here (core_gap, core_windings, scroll1_switch) is a nearest-point
distance on a grid whose nodes are ~900 um apart, so it carries a discretisation
floor of unknown size - dividing by it would launder that into a clean-looking
percentage. The coarse voxel is a denominator that is actually known.
"""
import json, os, sys, numpy as np, l0, xres

CACHE, FINE = xres.CACHE, f"{xres.CACHE}/_fine"
UM = xres.UM


def pairs_available():
    """Segments with both a coarse mesh and a 7.91 um one in the cache."""
    segs = {d.split("__")[0] for d in os.listdir(FINE) if "7.91um" in d}
    segs |= {d.split("-on-")[0] for d in os.listdir(FINE) if "7.91um" in d and "__" not in d}
    out = []
    for seg in sorted(segs):
        a, b = xres.find_mesh(seg, "45.5"), xres.find_mesh(seg, "7.91")
        if a and b: out.append((seg, a, b))
    return out


def one(seg, p45, p79, n=1500):
    P45, v45 = xres.load_fixed(p45, "45.5")
    P79, v79 = xres.load_fixed(p79, "7.91")
    Q = P45[v45]
    if len(Q) < 50: return None
    Q = Q[::max(1, len(Q) // n)]
    sd = xres.surface_dist(Q, P79, v79)
    ok = ~np.isnan(sd)
    Q, sd = Q[ok], sd[ok] * UM
    if len(sd) < 50: return None
    blk = np.floor(Q / 500).astype(np.int64)    # ~4 mm cells
    _, inv = np.unique(blk, axis=0, return_inverse=True)
    cnt = np.bincount(inv)
    big = [k for k in range(len(cnt)) if cnt[k] >= 15]
    mu = np.array([abs(sd[inv == k].mean()) for k in big]) if big else np.array([np.nan])
    sg = np.array([sd[inv == k].std() for k in big]) if big else np.array([np.nan])
    a = np.abs(sd)
    return dict(seg=seg, n=len(sd), grid45=list(P45.shape[:2]), grid79=list(P79.shape[:2]),
                centroid=[round(float(x), 1) for x in Q.mean(0)],
                p50=round(float(np.median(a)), 1), p90=round(float(np.percentile(a, 90)), 1),
                p99=round(float(np.percentile(a, 99)), 1), max=round(float(a.max()), 1),
                p50_coarse_vx=round(float(np.median(a) / 45.532), 2),
                block_mean_p50=round(float(np.median(mu)), 1),
                block_sd_p50=round(float(np.median(sg)), 1), blocks=len(big))


if __name__ == "__main__":
    todo = pairs_available()
    print(f"{len(todo)} segments with both a 45.5 um and a 7.91 um mesh", flush=True)
    res = []
    for i, (seg, a, b) in enumerate(todo):
        try:
            r = one(seg, a, b)
        except Exception as e:
            print(f"[{i+1}/{len(todo)}] {seg}: ERR {type(e).__name__}: {e}", flush=True); continue
        if r is None:
            print(f"[{i+1}/{len(todo)}] {seg}: too few points", flush=True); continue
        res.append(r)
        print(f"[{i+1}/{len(todo)}] {seg:<24} p50={r['p50']:6.1f} p90={r['p90']:6.1f} "
              f"p99={r['p99']:6.1f} um ({r['p50_coarse_vx']:.2f} coarse vx)  "
              f"blk_mean={r['block_mean_p50']:6.1f} blk_sd={r['block_sd_p50']:6.1f}", flush=True)
    json.dump(res, open("xres_results.json", "w"), indent=1)
    a = np.array([r["p50"] for r in res])
    print(f"\nwrote xres_results.json  ({len(res)} segments)")
    print(f"p50 of per-segment p50: {np.median(a):.1f} um = {np.median(a)/45.532:.2f} coarse voxels"
          f"   range {a.min():.1f}-{a.max():.1f} um")
