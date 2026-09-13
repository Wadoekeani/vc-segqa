"""List every Scroll 1 segment's coarse tifxyz and download it.

Lists each segment's mesh/ folder separately - listing the whole segments/
prefix would enumerate every surface-volume zarr chunk (hundreds of thousands).
"""
import os, sys, urllib.request, urllib.parse, xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"

def _get(url):
    return urllib.request.urlopen(url, timeout=60).read()

def list_dir(prefix):
    """One level: (subdir prefixes, keys)."""
    url = f"{BUCKET}/?list-type=2&delimiter=/&max-keys=1000&prefix={urllib.parse.quote(prefix)}"
    root = ET.fromstring(_get(url))
    return ([p.findtext(NS+"Prefix") for p in root.findall(NS+"CommonPrefixes")],
            [c.findtext(NS+"Key") for c in root.findall(NS+"Contents")])

def find_mesh(seg_prefix):
    """Pick the coarse (45um) tifxyz for one segment; skip intermediate/."""
    dirs, _ = list_dir(seg_prefix + "mesh/")
    cand = [d for d in dirs if d.rstrip("/").endswith(".tifxyz")]
    if not cand: return None
    coarse = [d for d in cand if "45.5" in d]
    return (coarse or cand)[0]

def fetch(args):
    d, out = args
    seg = d.split("/")[2]
    dest = os.path.join(out, seg)
    os.makedirs(dest, exist_ok=True)
    n = 0
    for f in ("x.tif", "y.tif", "z.tif", "meta.json"):
        p = os.path.join(dest, f)
        if os.path.exists(p) and os.path.getsize(p) > 0: continue
        try:
            urllib.request.urlretrieve(f"{BUCKET}/{d}{f}", p); n += 1
        except Exception as e:
            print(f"  {seg}/{f}: {e}")
    return seg, d, n

if __name__ == "__main__":
    out = sys.argv[1]
    segs, _ = list_dir("PHercParis4/segments/")
    print(f"{len(segs)} segments")
    with ThreadPoolExecutor(8) as ex:
        meshes = [(s, m) for s, m in zip(segs, ex.map(find_mesh, segs)) if m]
    print(f"{len(meshes)} have a tifxyz mesh "
          f"({sum('45.5' in m for _, m in meshes)} coarse)")
    with ThreadPoolExecutor(8) as ex:
        got = list(ex.map(fetch, [(m, out) for _, m in meshes]))
    print(f"downloaded {sum(n for _, _, n in got)} files into {out}")
    for s, _ in [(s, m) for s, m in zip(segs, [None]*len(segs))][:0]: pass
    missing = [s.split("/")[2] for s in segs if s not in [x[0] for x in meshes]]
    if missing: print(f"no mesh: {len(missing)} -> {missing[:5]}")
