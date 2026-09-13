"""Export per-vertex QA flags so a browser can colour the mesh without redoing
the volume sampling. One byte per grid cell: 0 inside, 1 outside, 2 no data.

Flags are subsampled to at most MAXG x MAXG - a 2.4 um mesh has millions of
cells, and the viewer only needs enough to show the shape of a failure. The
viewer applies the same stride when it parses the TIFF, so positions still
stream from S3 and nothing but these small files is served locally.
"""
import os, re, sys, json, numpy as np
from l0 import read_tifxyz_plane
from scan import seg_mesh_dirs, sample_mask, fetch_mesh, CACHE
from survey import BUCKET
from concurrent.futures import ThreadPoolExecutor

MAXG = 128

def stride_for(w, h):
    return max(1, int(np.ceil(max(w, h) / MAXG)))

def grid(d):
    P = np.dstack([read_tifxyz_plane(f"{d}/{a}.tif") for a in "xyz"])
    v = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P, v

def export(scroll, plan, root, out):
    vol, level = plan["volume"], plan["level"]
    dirs = seg_mesh_dirs(scroll, vol.split("-")[0])
    os.makedirs(out, exist_ok=True)
    with ThreadPoolExecutor(8) as ex:
        list(ex.map(fetch_mesh, [(d, f"{root}/{scroll}/{n}") for n, d in dirs.items()]))
    cache, index = {}, []
    for name in sorted(dirs):
        d = f"{root}/{scroll}/{name}"
        if not os.path.exists(f"{d}/x.tif"): continue
        try: P, v = grid(d)
        except Exception as e: print(f"   skip {name}: {e}"); continue
        h, w = v.shape
        flags = np.full((h, w), 2, np.uint8)
        pts = P[v]
        if len(pts):
            flags[v] = (sample_mask(scroll, vol, level, plan["shape"], pts, cache) == 0)
        s = stride_for(w, h)
        sub = np.ascontiguousarray(flags[::s, ::s])
        sub.tofile(f"{out}/{name}.bin")
        wn = re.search(r"w(\d{3})", name)
        index.append(dict(segment=name, w=w, h=h, stride=s,
                          sw=int(sub.shape[1]), sh=int(sub.shape[0]),
                          mesh=f"{BUCKET}/{dirs[name]}",
                          winding=int(wn.group(1)) if wn else None,
                          verts=int(v.sum()),
                          outside=float((flags[v] == 1).mean()) if v.any() else 0.0))
        print(f"   {name[:44]:44s} {w:>5}x{h:<5} /{s} {index[-1]['outside']:>7.2%}")
    index.sort(key=lambda r: (r["winding"] is None, r["winding"] or 0))
    tot = sum(r["verts"] for r in index) or 1
    meta = dict(scroll=scroll, volume=vol, level=level, um=plan["eff_um"],
                segments=index,
                outside=sum(r["verts"]*r["outside"] for r in index)/tot,
                verts=tot)
    json.dump(meta, open(f"{out}/index.json", "w"), indent=1)
    return meta

if __name__ == "__main__":
    plans = json.load(open(sys.argv[1]))
    want = sys.argv[2:] or [p["scroll"] for p in plans]
    root, out = CACHE, "viewer/data"
    summary = []
    for p in plans:
        if p["scroll"] not in want: continue
        print(f"=== {p['scroll']}")
        m = export(p["scroll"], p, root, f"{out}/{p['scroll']}")
        mb = sum(os.path.getsize(f"{out}/{p['scroll']}/{r['segment']}.bin")
                 for r in m["segments"]) / 1e6
        summary.append(dict(scroll=p["scroll"], segments=len(m["segments"]),
                            verts=m["verts"], outside=m["outside"], mb=round(mb, 2)))
        # rewrite after every scroll so the viewer works on partial output
        json.dump(summary, open(f"{out}/scrolls.json", "w"), indent=1)
        print(f"   -- {len(m['segments'])} segments, {mb:.2f} MB of flags")
    json.dump(summary, open(f"{out}/scrolls.json", "w"), indent=1)
    print(f"\n{len(summary)} scrolls -> {out}/scrolls.json "
          f"({sum(s['mb'] for s in summary):.1f} MB total)")
