"""Corpus-wide QA scan: do traced vertices land inside the scroll?

For every scroll with published segments and a masked volume, sample the mask
at each mesh vertex. Voxels outside the scroll body are exactly 0 in the
published masked volumes, so this needs no labels - the only quality signal
available without ground truth.
"""
import os, re, sys, json, urllib.request, urllib.error, numpy as np
from concurrent.futures import ThreadPoolExecutor
from survey import list_dir, BUCKET
from l0 import read_tifxyz_plane

CH = 128
CACHE = os.environ.get("SEGQA_CACHE", os.path.expanduser("~/.cache/segqa"))

def seg_mesh_dirs(scroll, volid):
    """Every segment's tifxyz registered to this volume."""
    segs, _ = list_dir(f"{scroll}/segments/")
    def one(sp):
        try: dirs, _ = list_dir(sp + "mesh/")
        except Exception: return None
        for d in dirs:
            if re.search(rf"-on-{volid}-[\d.]+um\.tifxyz/$", d): return d
        return None
    with ThreadPoolExecutor(16) as ex:
        return {sp.rstrip("/").split("/")[-1]: d
                for sp, d in zip(segs, ex.map(one, segs)) if d}

def fetch_mesh(args):
    d, dest = args
    os.makedirs(dest, exist_ok=True)
    for a in "xyz":
        p = f"{dest}/{a}.tif"
        if os.path.exists(p) and os.path.getsize(p): continue
        try: urllib.request.urlretrieve(f"{BUCKET}/{d}{a}.tif", p)
        except Exception as e: return f"{dest}: {e}"
    return None

def load_pts(d):
    P = np.dstack([read_tifxyz_plane(f"{d}/{a}.tif") for a in "xyz"])
    v = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P[v]

def get_chunk(url):
    try:
        raw = urllib.request.urlopen(url, timeout=120).read()
        return np.frombuffer(raw, np.uint8).reshape(CH, CH, CH)
    except Exception:
        return None                      # missing chunk == all fill_value == outside

def sample_mask(scroll, vol, level, shape, pts, cache):
    """Fetch only the chunks the vertices actually touch."""
    idx = np.rint(pts[:, ::-1] / (1 << level)).astype(np.int64)     # xyz -> zyx
    inb = np.all((idx >= 0) & (idx < np.array(shape)), axis=1)
    out = np.zeros(len(pts), np.uint8)
    if not inb.any(): return out
    ii = idx[inb]
    ci = ii // CH
    # one integer per chunk so grouping is a single np.unique, not a loop per chunk
    nb = [(s + CH - 1)//CH for s in shape]
    code = (ci[:, 0]*nb[1] + ci[:, 1])*nb[2] + ci[:, 2]
    keys, inv = np.unique(code, return_inverse=True)
    base = f"{BUCKET}/{scroll}/volumes/{vol}/{level}"
    coords = [(int(k)//(nb[1]*nb[2]), int(k)//nb[2] % nb[1], int(k) % nb[2]) for k in keys]
    todo = [c for c in coords if c not in cache]
    if todo:
        with ThreadPoolExecutor(12) as ex:
            for c, blk in zip(todo, ex.map(get_chunk, [f"{base}/{z}/{y}/{x}" for z, y, x in todo])):
                cache[c] = blk
    loc = ii % CH
    vals = np.zeros(len(ii), np.uint8)
    for gi, c in enumerate(coords):
        blk = cache[c]
        if blk is None: continue
        m = inv == gi
        l = loc[m]
        vals[m] = blk[l[:, 0], l[:, 1], l[:, 2]]
    out[inb] = vals
    return out

def scan(plan, root, results):
    scroll, vol, level = plan["scroll"], plan["volume"], plan["level"]
    volid = vol.split("-")[0]
    dirs = seg_mesh_dirs(scroll, volid)
    print(f"\n=== {scroll}: {len(dirs)}/{plan['segments']} segments on {vol} "
          f"(L{level}, {plan['eff_um']} um)")
    out = f"{root}/{scroll}"
    with ThreadPoolExecutor(8) as ex:
        for err in ex.map(fetch_mesh, [(d, f"{out}/{s}") for s, d in dirs.items()]):
            if err: print("  ", err)
    cache, rows = {}, []
    for s in sorted(dirs):
        try: pts = load_pts(f"{out}/{s}")
        except Exception as e: print(f"   skip {s}: {e}"); continue
        if not len(pts): continue
        v = sample_mask(scroll, vol, level, plan["shape"], pts, cache)
        rows.append(dict(scroll=scroll, segment=s, verts=len(pts),
                         outside=float((v == 0).mean()),
                         median=float(np.median(v[v > 0])) if (v > 0).any() else 0.0))
    rows.sort(key=lambda r: -r["outside"])
    for r in rows[:5]:
        print(f"   {r['segment'][:44]:44s} {r['verts']:>8,} {r['outside']:>7.2%}")
    tot = sum(r["verts"] for r in rows)
    bad = sum(r["verts"]*r["outside"] for r in rows)
    print(f"   -- {len(rows)} segments, {tot:,} verts, {bad/max(tot,1):.2%} outside "
          f"({sum(r['outside'] > .05 for r in rows)} segments over 5%)")
    results += rows

if __name__ == "__main__":
    plans = json.load(open(sys.argv[1]))
    root = sys.argv[2] if len(sys.argv) > 2 else CACHE
    results = []
    for p in plans:
        try: scan(p, root, results)
        except Exception as e: print(f"{p['scroll']}: FAILED {e}")
    json.dump(results, open("scan_results.json", "w"), indent=1)
    print(f"\nwrote scan_results.json: {len(results)} segments across "
          f"{len({r['scroll'] for r in results})} scrolls")
