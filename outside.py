"""Do traced vertices actually land inside the scroll?

The published volume is masked: voxels outside the scroll body are exactly 0.
That makes "is this vertex on the object at all" an objective check needing no
labels - the one QA signal available without ground truth.
"""
import os, sys, glob, urllib.request, numpy as np
from concurrent.futures import ThreadPoolExecutor
from coverage import load_seg

BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
VOL = "PHercParis4/volumes/20260310170716-45.532um-11.0m-74keV-masked.zarr"
CH, LEVEL, SHAPE = 128, 3, (509, 283, 283)          # 364 um/voxel - plenty for a mask

def get_chunk(t):
    zc, yc, xc = t
    try:
        raw = urllib.request.urlopen(f"{BUCKET}/{VOL}/{LEVEL}/{zc}/{yc}/{xc}", timeout=120).read()
        return t, np.frombuffer(raw, np.uint8).reshape(CH, CH, CH)
    except Exception:
        return t, np.zeros((CH, CH, CH), np.uint8)

def load_volume(cache):
    if os.path.exists(cache):
        return np.load(cache)["v"]
    nz, ny, nx = [(s + CH - 1)//CH for s in SHAPE]
    jobs = [(z, y, x) for z in range(nz) for y in range(ny) for x in range(nx)]
    vol = np.zeros((nz*CH, ny*CH, nx*CH), np.uint8)
    with ThreadPoolExecutor(16) as ex:
        for (zc, yc, xc), c in ex.map(get_chunk, jobs):
            vol[zc*CH:(zc+1)*CH, yc*CH:(yc+1)*CH, xc*CH:(xc+1)*CH] = c
    vol = vol[:SHAPE[0], :SHAPE[1], :SHAPE[2]]
    np.savez_compressed(cache, v=vol)
    return vol

def sample(vol, pts):
    """Nearest-neighbour sample; points off the array count as outside."""
    idx = np.rint(pts[:, ::-1] / (1 << LEVEL)).astype(int)      # xyz -> zyx
    ok = np.all((idx >= 0) & (idx < np.array(SHAPE)), axis=1)
    out = np.zeros(len(pts), np.uint8)
    out[ok] = vol[idx[ok, 0], idx[ok, 1], idx[ok, 2]]
    return out, ok

if __name__ == "__main__":
    root, cache = sys.argv[1], sys.argv[2]
    vol = load_volume(cache)
    print(f"volume {vol.shape} level {LEVEL}, inside mask {np.count_nonzero(vol)/vol.size:.1%}")
    segs = {os.path.basename(d.rstrip("/")): load_seg(d)
            for d in sorted(glob.glob(f"{root}/*/"))}
    rows = []
    for n, p in segs.items():
        v, ok = sample(vol, p)
        rows.append((float((v == 0).mean()), n, len(p), float(np.median(v[v > 0])) if (v > 0).any() else 0))
    rows.sort(reverse=True)
    print(f"\n{'segment':40s} {'verts':>8} {'outside mask':>13} {'median val':>11}")
    for frac, n, np_, med in rows:
        flag = "  <-- " + ("BAD" if frac > .05 else "check") if frac > .01 else ""
        print(f"{n:40s} {np_:>8,} {frac:>12.2%} {med:>11.0f}{flag}")
    tot = np.concatenate([sample(vol, p)[0] for p in segs.values()])
    print(f"\nall segments: {(tot==0).mean():.2%} of {len(tot):,} vertices land outside the scroll")
