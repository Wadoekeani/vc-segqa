"""L0 metrics on a tifxyz quad grid. Pure numpy, no deps."""
import struct, numpy as np

def _lzw(data, expect):
    """TIFF-flavoured LZW: MSB-first, and the code width grows one code earlier
    than plain LZW ("early change") - the usual way a hand-rolled decoder ends
    up subtly wrong."""
    out = bytearray()
    table = [bytes([i]) for i in range(256)] + [b"", b""]
    prev, width, pos, nbits = None, 9, 0, len(data) * 8
    while pos + width <= nbits and len(out) < expect:
        byte, bit = divmod(pos, 8)
        chunk = int.from_bytes(data[byte:byte + 3].ljust(3, b"\0"), "big")
        code = (chunk >> (24 - bit - width)) & ((1 << width) - 1)
        pos += width
        if code == 256:
            del table[258:]; width, prev = 9, None; continue
        if code == 257: break
        if code < len(table) and (code < 256 or table[code]):
            entry = table[code]
        elif prev is not None:
            entry = prev + prev[:1]
        else:
            raise ValueError("corrupt LZW stream")
        out += entry
        if prev is not None: table.append(prev + entry[:1])
        prev = entry
        if len(table) + 1 >= (1 << width) and width < 12: width += 1
    return bytes(out)


def _undo_predictor(buf, tw, th, predictor, dt):
    """Reverse TIFF prediction over a tile of tw x th float32 samples."""
    if predictor in (None, 1):
        return np.frombuffer(buf, dt, tw * th).reshape(th, tw)
    a = np.frombuffer(buf, np.uint8, tw * th * 4).reshape(th, tw * 4).astype(np.int64)
    a = np.cumsum(a, axis=1).astype(np.uint8)          # horizontal differencing
    if predictor == 2:
        return a.view(dt).reshape(th, tw)
    if predictor == 3:
        # bytes are stored plane by plane along the row, most significant first
        planes = a.reshape(th, 4, tw)
        return np.ascontiguousarray(planes.transpose(0, 2, 1)).view(">f4").reshape(th, tw)
    raise ValueError(f"unsupported predictor {predictor}")


def read_tifxyz_plane(path):
    """One plane of a tifxyz mesh.

    Most published meshes are classic TIFF, uncompressed, single strip - the
    whole reader for those is two lines. A few are BigTIFF, tiled, LZW with the
    floating-point predictor, so those paths exist too; a classic-only reader
    fails on them with a struct error that names nothing.
    """
    d = open(path, 'rb').read()
    bo = '<' if d[:2] == b'II' else '>'
    magic, = struct.unpack(bo + 'H', d[2:4])
    if magic == 42:
        off, = struct.unpack(bo + 'I', d[4:8])
        n, = struct.unpack(bo + 'H', d[off:off + 2])
        ent, esz, vo = off + 2, 12, 8
        widths = {1: 'B', 3: 'H', 4: 'I'}
    elif magic == 43:
        offsz, zero = struct.unpack(bo + 'HH', d[4:8])
        assert offsz == 8 and zero == 0, f"odd BigTIFF header in {path}"
        off, = struct.unpack(bo + 'Q', d[8:16])
        n, = struct.unpack(bo + 'Q', d[off:off + 8])
        ent, esz, vo = off + 8, 8 + 12, 12
        widths = {1: 'B', 3: 'H', 4: 'I', 16: 'Q', 17: 'q'}
    else:
        raise ValueError(f"not a TIFF (magic {magic}): {path}")

    t = {}
    for i in range(n):
        e = ent + i * esz
        tag, typ = struct.unpack(bo + 'HH', d[e:e + 4])
        w = widths.get(typ)
        if w: t[tag] = struct.unpack(bo + w, d[e + vo:e + vo + struct.calcsize(w)])[0]

    w, h = t[256], t[257]
    assert t[258] == 32 and t[339] == 3, f"expect float32, got bits={t[258]} fmt={t[339]}"
    dt = '<f4' if bo == '<' else '>f4'
    comp, pred = t[259], t.get(317)

    if 324 in t:                                   # tiled
        tw, th = t[322], t[323]
        raw = d[t[324]:t[324] + t[325]]
        if comp == 5: raw = _lzw(raw, tw * th * 4)
        elif comp != 1: raise ValueError(f"unsupported compression {comp} in {path}")
        tile = _undo_predictor(raw, tw, th, pred, dt)
        assert tw >= w and th >= h, f"{path}: multi-tile images not handled"
        return np.ascontiguousarray(tile[:h, :w]).astype('<f4')

    raw = d[t[273]:t[273] + t[279]]
    if comp == 5: raw = _lzw(raw, w * h * 4)
    elif comp != 1: raise ValueError(f"unsupported compression {comp} in {path}")
    return _undo_predictor(raw, w, h, pred, dt).astype('<f4')


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

def _tiff_selfcheck():
    """Decode a stream small enough to verify by hand.

    Codes (9 bits each): CLEAR, 'A', 'B', 258, EOI. After emitting A then B the
    decoder must have learned 258 = "AB", so code 258 emits "AB" and the whole
    stream is "ABAB". This exercises the table growth and the deferred entry,
    which is where a hand-written LZW usually goes wrong.
    """
    codes = [256, ord("A"), ord("B"), 258, 257]
    acc = nacc = 0
    buf = bytearray()
    for c in codes:                       # plain 9-bit MSB-first packing
        acc = (acc << 9) | c; nacc += 9
        while nacc >= 8:
            nacc -= 8; buf.append((acc >> nacc) & 0xFF)
    if nacc: buf.append((acc << (8 - nacc)) & 0xFF)
    got = _lzw(bytes(buf), 16)
    assert got == b"ABAB", f"LZW decoded {got!r}, expected b'ABAB'"

    # predictor 2 is plain horizontal differencing over bytes
    a = np.arange(12, dtype=np.uint8).reshape(1, 12)
    diff = np.diff(a, axis=1, prepend=0).astype(np.uint8).tobytes()
    assert _undo_predictor(diff, 3, 1, 2, "<f4").tobytes() == a.tobytes(), "predictor 2 failed"

    # predictor 3 regroups the bytes of each float into planes, MSB first
    vals = np.array([[1.5, -2.25, 1e4]], dtype=">f4")
    planes = vals.view(np.uint8).reshape(1, 3, 4).transpose(0, 2, 1).reshape(1, 12)
    enc = np.diff(planes.astype(np.int64), axis=1, prepend=0).astype(np.uint8).tobytes()
    out = _undo_predictor(enc, 3, 1, 3, ">f4")
    assert np.allclose(out, vals), f"predictor 3 gave {out}"
    print("tiff self-check ok (LZW table growth, predictor 2 and 3)")


if __name__ == "__main__":
    import sys
    _tiff_selfcheck()
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
