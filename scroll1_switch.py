"""Adjacent-winding check on Scroll 1, using windings recovered by unroll.py.

Same comparison that found PHerc0139 w46 - the range-named segments are first
split into one point set per winding, so real adjacent pairs exist.
"""
import os, sys, glob, json, numpy as np
from collections import defaultdict
from unroll import slices, range_of
from neighbour import nn_dist
from coverage import load_seg
from axis import fit_axis
from scan import CACHE

UM = 45.532

def build(scroll, batch):
    dirs = [d for d in sorted(glob.glob(f"{CACHE}/{scroll}/*/"))
            if os.path.exists(f"{d}/x.tif")
            and range_of(os.path.basename(d.rstrip("/")))
            and os.path.basename(d.rstrip("/")).startswith(batch)]
    allp = np.concatenate([load_seg(d) for d in dirs])
    zc, cen = fit_axis(allp)
    by = defaultdict(list)
    for d in dirs:
        sl, msg = slices(d, zc, cen)
        if not sl:
            print(f"  skip {os.path.basename(d.rstrip('/'))}: {msg}"); continue
        for k, pts in sl.items(): by[k].append(pts)
    return {k: np.concatenate(v) for k, v in by.items()}

def run(by, label):
    ws = sorted(by)
    pairs = [(a, a + 1) for a in ws if a + 1 in by]
    print(f"\n=== {label}: windings {min(ws)}–{max(ws)}, {len(pairs)} adjacent pairs")
    rng = np.random.default_rng(0)
    rows, gaps = [], []
    for a, b in pairs:
        A, B = by[a], by[b]
        sub = A[rng.choice(len(A), min(1800, len(A)), replace=False)]
        d = nn_dist(sub, B, cell=max(8.0, 400.0 / UM))
        ov = float(np.isfinite(d).mean()); d = d[np.isfinite(d)]
        if len(d) < 50 or ov < 0.5:
            print(f"  {a:>4}/{b:<4} skipped: {ov:.0%} overlap"); continue
        q = np.percentile(d, [10, 50, 90])
        rows.append(dict(a=a, b=b, overlap=ov, p50=float(q[1]), p10=float(q[0]),
                         p90=float(q[2]), um=float(q[1] * UM)))
        gaps.append(q[1])
    if not rows: return []
    med = float(np.median(gaps))
    print(f"  typical adjacent-winding gap: {med:.1f} vx = {med*UM:.0f} um")
    print(f"  {'pair':>9} {'ovl':>5} {'p50':>7} {'um':>6} {'p50/med':>8}")
    for r in rows:
        r["ratio"] = r["p50"] / med
        r["flag"] = ("COINCIDENT" if r["ratio"] < 0.35 else
                     "GAP" if r["ratio"] > 1.8 else "")
        print(f"  {r['a']:>4}/{r['b']:<4} {r['overlap']:>5.0%} {r['p50']:>7.1f} "
              f"{r['um']:>6.0f} {r['ratio']:>8.2f}"
              + (f"  <-- {r['flag']}" if r["flag"] else ""))
    return rows

if __name__ == "__main__":
    out = {}
    for batch in ("20260623", "20260701"):
        by = build("PHercParis4", batch)
        if by: out[batch] = run(by, f"PHercParis4 batch {batch}")
    json.dump(out, open("scroll1_switch.json", "w"), indent=1)
    tot = sum(len(v) for v in out.values())
    bad = sum(1 for v in out.values() for r in v if r["flag"])
    print(f"\nwrote scroll1_switch.json: {tot} pairs, {bad} flagged")
