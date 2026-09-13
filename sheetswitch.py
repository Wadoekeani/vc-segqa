"""Run the adjacent-winding check over a scroll and flag suspicious pairs."""
import os, re, sys, json, glob, numpy as np
from coverage import load_seg
from neighbour import nn_dist
from scan import CACHE

def winding(name):
    m = re.findall(r"w(\d{3})", name)
    return int(m[0]) if len(m) == 1 else None

def segments(scroll):
    out = {}
    for d in sorted(glob.glob(f"{CACHE}/{scroll}/*/")):
        n = os.path.basename(d.rstrip("/"))
        w = winding(n)
        if w is None or not os.path.exists(f"{d}/x.tif"): continue
        out.setdefault(w, []).append(d)
    return out

def run(scroll, um, nsample=1800):
    segs = segments(scroll)
    ws = sorted(segs)
    pairs = [(a, a + 1) for a in ws if a + 1 in segs]
    if not pairs: return []
    print(f"\n=== {scroll}  {len(pairs)} adjacent pairs  ({um} um/voxel)")
    # 2.4 um meshes run to millions of vertices each; cap both sides and keep
    # only the two windings in play so memory stays flat across a long run.
    MAXB = 600_000
    pts, rng0 = {}, np.random.default_rng(1)
    def get(w):
        if w not in pts:
            P = np.concatenate([load_seg(d) for d in segs[w]])
            if len(P) > MAXB:
                P = P[rng0.choice(len(P), MAXB, replace=False)]
            pts[w] = P
        for k in list(pts):
            if k < w - 1: del pts[k]
        return pts[w]

    rows, gaps = [], []
    rng = np.random.default_rng(0)
    for a, b in pairs:
        A, B = get(a), get(b)
        if len(A) < 50 or len(B) < 50: continue
        idx = rng.choice(len(A), min(nsample, len(A)), replace=False)
        sub = A[idx]
        cell = max(8.0, 400.0 / um)
        # cell a bit over the expected spacing so the 27-cell search always covers it
        d = nn_dist(sub, B, cell=cell)
        overlap = float(np.isfinite(d).mean())
        d = d[np.isfinite(d)]
        if len(d) < 50 or overlap < 0.5: 
            print(f"  {a:>4}/{b:<4} skipped: only {overlap:.0%} of sampled points "
                  f"have a neighbour (segments barely overlap)")
            continue
        # Sampling density puts a floor under any nearest-neighbour distance, and
        # the meshes differ in density by 10x across scrolls. Measure it on the same
        # points against their own mesh and report it, so a pair whose gap is not
        # actually resolvable can be recognised rather than trusted. It is reported,
        # not thresholded: the tifxyz grid step (~20 vx) is wider than the sheet
        # spacing (~16 vx), so gap < floor is the normal, healthy case.
        rest = np.setdiff1d(np.arange(len(A)), idx, assume_unique=False)
        if len(rest) > len(B): rest = rng.choice(rest, len(B), replace=False)
        floor = float(np.median(nn_dist(sub, A[rest], cell=cell)))
        q = np.percentile(d, [10, 50, 90])
        rows.append(dict(scroll=scroll, a=a, b=b, n=len(d), overlap=overlap,
                         p10=float(q[0]), p50=float(q[1]), p90=float(q[2]),
                         floor=float(floor), snr=float(q[1] / max(floor, 1e-9)),
                         um=float(q[1] * um)))
        gaps.append(q[1])
    med = float(np.median(gaps))
    print(f"  typical adjacent-winding gap: {med:.1f} vx = {med*um:.0f} um")
    print(f"  {'pair':>9} {'ovl':>5} {'floor':>7} {'p50':>7} {'gap/floor':>10} {'p50/med':>8}")
    for r in rows:
        ratio = r["p50"] / med
        flag = ""
        if ratio < 0.35:  flag = "  <-- COINCIDENT (same sheet?)"
        elif ratio > 1.8: flag = "  <-- GAP (winding skipped?)"
        r["ratio"] = ratio; r["flag"] = flag.strip(" <-")
        print(f"  {r['a']:>4}/{r['b']:<4} {r['overlap']:>5.0%} {r['floor']:>7.1f} "
              f"{r['p50']:>7.1f} {r['snr']:>10.2f} {ratio:>8.2f}{flag}")
    return rows

if __name__ == "__main__":
    plans = {p["scroll"]: p for p in json.load(open("scan_plan.json"))}
    out = []
    for s in sys.argv[1:]:
        out += run(s, plans[s]["base_um"])
    json.dump(out, open("sheetswitch.json", "w"), indent=1)
    print(f"\nwrote sheetswitch.json ({len(out)} pairs)")
