"""Fetch what xres_scan.py needs: the shared-frame transforms and the meshes cut
on a second volume.

fetch_meshes.py already pulls the coarse (45.5 um) mesh of every segment. This
pulls the finer ones plus each volume's transform.json, and names the cache dir
`<segment>__<mesh dir name>` because a mesh folder's prefix is not always the
segment id (`20260603222816-20231005123336_v2_flatboi` ships
`20260603222816-on-...`), so the segment has to be kept alongside.
"""
import os, sys, urllib.request, urllib.parse, xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from fetch_meshes import BUCKET, list_dir

SCROLL = "PHercParis4"
CACHE = os.path.expanduser(os.environ.get("SEGQA_CACHE", "~/.cache/segqa"))
VOL_ZARRS = ["20260310170716-45.532um-11.0m-74keV-masked",
             "20260411134726-2.400um-0.2m-78keV-masked",
             "20260608103018-1.129um-0.2m-78keV-masked"]


def fetch_transforms():
    d = f"{CACHE}/_transforms"; os.makedirs(d, exist_ok=True)
    for v in VOL_ZARRS:
        p = f"{d}/{v}.json"
        if os.path.exists(p) and os.path.getsize(p): continue
        try:
            urllib.request.urlretrieve(f"{BUCKET}/{SCROLL}/volumes/{v}.zarr/transform.json", p)
            print(f"  transform {v[:14]}")
        except Exception as e:
            print(f"  transform {v[:14]}: {e}")     # the 110keV 45.5 um volume has none


def mesh_dirs(tag):
    segs, _ = list_dir(f"{SCROLL}/segments/")
    def one(p):
        try: return [d for d in list_dir(p + "mesh/")[0] if tag in d and d.rstrip("/").endswith(".tifxyz")]
        except Exception: return []
    with ThreadPoolExecutor(12) as ex:
        return [d for r in ex.map(one, segs) for d in r]


def fetch_mesh(d):
    seg, name = d.split("/")[2], d.rstrip("/").split("/")[-1]
    dest = f"{CACHE}/_fine/{seg}__{name}"
    tmp = dest + ".part"
    if all(os.path.exists(f"{dest}/{f}") and os.path.getsize(f"{dest}/{f}")
           for f in ("x.tif", "y.tif", "z.tif", "meta.json")):
        return None
    os.makedirs(tmp, exist_ok=True)
    try:
        for f in ("x.tif", "y.tif", "z.tif", "meta.json"):
            urllib.request.urlretrieve(f"{BUCKET}/{d}{f}", f"{tmp}/{f}")
    except Exception as e:
        return f"{seg}: {e}"
    os.replace(tmp, dest)                 # only a complete dir gets the real name
    return None


if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "7.91um"
    print("transforms:")
    fetch_transforms()
    dirs = mesh_dirs(tag)
    print(f"{len(dirs)} meshes matching {tag!r}")
    with ThreadPoolExecutor(8) as ex:
        errs = [e for e in ex.map(fetch_mesh, dirs) if e]
    for e in errs: print("  ERR", e)
    print(f"done ({len(dirs) - len(errs)} present, {len(errs)} failed) in {CACHE}/_fine")
