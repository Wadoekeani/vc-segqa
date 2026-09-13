"""How much does the escape rate depend on the sampling resolution?

The corpus scan samples each scroll at whatever pyramid level lands nearest
320 um, which differs between scrolls because the base voxel sizes do. If the
measurement drifts with level, cross-scroll numbers mean nothing. So measure
the drift instead of assuming it: same segments, several levels.
"""
import sys, json, urllib.request, numpy as np
from concurrent.futures import ThreadPoolExecutor
from scan import seg_mesh_dirs, fetch_mesh, load_pts, sample_mask, CACHE
from survey import BUCKET

def zshape(scroll, vol, level):
    try:
        za = json.load(urllib.request.urlopen(
            f"{BUCKET}/{scroll}/volumes/{vol}/{level}/.zarray", timeout=30))
        return za["shape"]
    except Exception:
        return None

if __name__ == "__main__":
    plan = next(p for p in json.load(open(sys.argv[1])) if p["scroll"] == sys.argv[2])
    levels = [int(a) for a in sys.argv[3].split(",")]
    nseg = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    scroll, vol, base = plan["scroll"], plan["volume"], plan["base_um"]

    dirs = seg_mesh_dirs(scroll, vol.split("-")[0])
    names = sorted(dirs)
    pick = [names[round(i * (len(names)-1) / (nseg-1))] for i in range(nseg)]
    with ThreadPoolExecutor(6) as ex:
        list(ex.map(fetch_mesh, [(dirs[n], f"{CACHE}/{scroll}/{n}") for n in pick]))
    pts = {n: load_pts(f"{CACHE}/{scroll}/{n}") for n in pick}

    print(f"{scroll}  base {base} um  ({len(pick)} segments)\n")
    hdr = "".join(f"{'L%d' % L:>9}" for L in levels)
    print(f"{'segment':40s}{hdr}")
    um = {}
    table = {}
    for L in levels:
        shape = zshape(scroll, vol, L)
        if shape is None: print(f"  (no level {L})"); continue
        um[L] = round(base * (1 << L), 1)
        cache = {}
        for n in pick:
            v = sample_mask(scroll, vol, L, shape, pts[n], cache)
            table.setdefault(n, {})[L] = float((v == 0).mean())
    for n in pick:
        row = "".join(f"{table[n].get(L, float('nan'))*100:>8.2f}%" for L in levels if L in um)
        print(f"{n[:40]:40s}{row}")
    print(f"{'voxel size (um)':40s}" + "".join(f"{um[L]:>9.1f}" for L in levels if L in um))

    ok = [L for L in levels if L in um]
    tot = {L: np.mean([table[n][L] for n in pick]) for L in ok}
    print(f"\n{'mean escape rate':40s}" + "".join(f"{tot[L]*100:>8.2f}%" for L in ok))
    lo, hi = min(tot.values()), max(tot.values())
    print(f"\nspread across {um[ok[0]]:.0f}-{um[ok[-1]]:.0f} um: "
          f"{lo*100:.2f}% .. {hi*100:.2f}%  (factor {hi/max(lo,1e-9):.2f}x)")
