"""L0 metrics on a tifxyz quad grid. Pure numpy, no deps."""
import struct, numpy as np

def read_tifxyz_plane(path):
    d = open(path, 'rb').read()
    bo = '<' if d[:2] == b'II' else '>'
    off, = struct.unpack(bo+'I', d[4:8])
    n, = struct.unpack(bo+'H', d[off:off+2])
    t = {}
    for i in range(n):
        e = off+2+i*12
        tag, typ, cnt = struct.unpack(bo+'HHI', d[e:e+8])
        t[tag] = struct.unpack(bo+'H' if typ == 3 else bo+'I', d[e+8:e+10 if typ == 3 else e+12])[0]
    w, h = t[256], t[257]
    assert t[258] == 32 and t[339] == 3 and t[259] == 1, "expect uncompressed float32"
    return np.frombuffer(d, '<f4', w*h, t[273]).reshape(h, w)

def load(dirpath):
    P = np.dstack([read_tifxyz_plane(f"{dirpath}/{a}.tif") for a in "xyz"])
    valid = np.isfinite(P).all(-1) & (P != -1).any(-1) & (P != 0).all(-1)
    return P, valid

def l0(P, valid):
    # edge lengths along the two grid directions
    du = np.linalg.norm(np.diff(P, axis=1), axis=-1)   # (h, w-1)
    dv = np.linalg.norm(np.diff(P, axis=0), axis=-1)   # (h-1, w)
    qu, qv = du[:-1], dv[:, :-1]                        # per-quad, aligned
    ok = valid[:-1, :-1] & valid[:-1, 1:] & valid[1:, :-1] & valid[1:, 1:]
    med = np.median(np.concatenate([qu[ok], qv[ok]]))
    # isotropy: 1.0 = perfect square quad; grid is uniform in param space
    aniso = np.where(ok, np.maximum(qu, qv) / np.maximum(np.minimum(qu, qv), 1e-9), np.nan)
    # stretch relative to the segment's own median edge
    stretch = np.where(ok, np.maximum(qu, qv) / med, np.nan)
    return dict(grid=P.shape[:2], valid_frac=valid.mean(), quads=int(ok.sum()),
                median_edge_vx=float(med), aniso=aniso, stretch=stretch)

def report(name, dirpath):
    P, valid = load(dirpath)
    m = l0(P, valid)
    a, s = m['aniso'][~np.isnan(m['aniso'])], m['stretch'][~np.isnan(m['stretch'])]
    print(f"\n== {name}  grid={m['grid']}  valid={m['valid_frac']:.1%}  quads={m['quads']}")
    print(f"   median edge {m['median_edge_vx']:.2f} vx")
    print(f"   anisotropy  p50={np.percentile(a,50):.2f} p95={np.percentile(a,95):.2f} p99.9={np.percentile(a,99.9):.2f} max={a.max():.1f}")
    print(f"   stretch     p50={np.percentile(s,50):.2f} p95={np.percentile(s,95):.2f} p99.9={np.percentile(s,99.9):.2f} max={s.max():.1f}")
    print(f"   quads >5x median stretch: {(s>5).sum()} ({(s>5).mean():.3%})")
    return m

if __name__ == "__main__":
    import sys
    # self-check first: runs with no data, fails loudly if the maths breaks
    u, v = np.meshgrid(np.arange(40.), np.arange(30.))
    flat = np.dstack([u, v, np.zeros_like(u)])
    c = l0(flat, np.ones(u.shape, bool))
    assert abs(c["median_edge_vx"] - 1.0) < 1e-6, c["median_edge_vx"]
    assert np.nanmax(c["aniso"]) < 1.0001, np.nanmax(c["aniso"])
    torn = flat.copy(); torn[15:, :, 1] += 50.0          # tear one row away
    assert np.nanmax(l0(torn, np.ones(u.shape, bool))["stretch"]) > 20
    print("self-check ok (flat grid -> 1.0, torn grid -> flagged)")
    for d in sys.argv[1:]:
        report(d, d)
