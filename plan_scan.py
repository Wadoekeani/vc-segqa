"""For each scroll: pick the coarsest masked volume that meshes are registered
to, and the pyramid level closest to a target voxel size."""
import re, sys, json, urllib.request, urllib.parse, numpy as np
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from survey import list_dir, BUCKET

TARGET_UM = 320.0

def zarray(vol, level):
    try:
        return json.load(urllib.request.urlopen(f"{BUCKET}/{vol}/{level}/.zarray", timeout=30))
    except Exception:
        return None

def pick(scroll):
    segs, _ = list_dir(f"{scroll}/segments/")
    vols, _ = list_dir(f"{scroll}/volumes/")
    masked = {v.rstrip("/").split("/")[-1]: v for v in vols if "masked" in v}
    def mesh_dirs(sp):
        try: return [d.rstrip("/").split("/")[-1] for d in list_dir(sp+"mesh/")[0]]
        except Exception: return []
    with ThreadPoolExecutor(16) as ex:
        allm = list(ex.map(mesh_dirs, segs))
    use = {}
    for ms in allm:
        for m in ms:
            g = re.search(r"-on-(\d+)-([\d.]+)um\.tifxyz$", m)
            if g: use.setdefault((g.group(1), float(g.group(2))), 0)
            if g: use[(g.group(1), float(g.group(2)))] += 1
    cand = [(res, vid, n) for (vid, res), n in use.items()
            if any(k.startswith(vid) for k in masked)]
    if not cand: return None
    res, vid, n = max(cand)                       # coarsest volume with meshes on it
    vname = next(k for k in masked if k.startswith(vid))
    levels = []
    for L in range(8):
        za = zarray(f"{scroll}/volumes/{vname}", L)
        if za is None: break
        levels.append((L, za["shape"], res * (1 << L)))
    best = min(levels, key=lambda t: abs(t[2] - TARGET_UM))
    mb = np.prod(best[1]) / 1e6
    return dict(scroll=scroll, segments=len(segs), meshes=n, volume=vname,
                base_um=res, level=best[0], shape=best[1],
                eff_um=round(best[2], 1), mb=round(mb, 1))

if __name__ == "__main__":
    out = []
    for s in sys.argv[1:]:
        try:
            r = pick(s)
        except Exception as e:
            print(f"{s}: {e}"); continue
        if not r: print(f"{s}: no masked volume with meshes"); continue
        out.append(r)
        print(f"{r['scroll']:14s} segs={r['segments']:>3} vol={r['base_um']:>7.3f}um "
              f"L{r['level']} -> {r['eff_um']:>6.1f}um  shape={r['shape']}  {r['mb']:>7.1f} MB")
    json.dump(out, open("scan_plan.json", "w"), indent=1)
    print(f"\ntotal segments {sum(r['segments'] for r in out)}, "
          f"volume download {sum(r['mb'] for r in out)/1000:.2f} GB")
