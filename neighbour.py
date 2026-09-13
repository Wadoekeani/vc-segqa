"""Sheet-switch detection by comparing segments of adjacent windings.

Neither of the other checks can see a sheet switch: the surface stays inside
the scroll (so the mask says nothing) and the parameterisation stays locally
isometric (so distortion says nothing). The signal only exists *between*
segments.

Segment names carry official winding numbers. Two segments one winding apart
must sit exactly one papyrus thickness apart, everywhere they overlap. So:

  distance collapses to ~0  -> both traced the same sheet; one of them jumped
  distance roughly doubles  -> a winding was skipped
  a step in the distance    -> the jump happened at that place

The winding labels supply the expectation, so this needs no quality labels
either - the same trick as the mask check.
"""
import numpy as np

def nn_dist(A, B, cell, want_index=False):
    """Distance from every point of A to the nearest point of B.

    A uniform spatial hash rather than a KD-tree: no scipy dependency, and the
    points are a thin surface so the grid stays sparse.
    """
    lo = B.min(0) - cell
    keyB = np.floor((B - lo) / cell).astype(np.int64)
    dims = keyB.max(0) + 2
    flat = (keyB[:, 0] * dims[1] + keyB[:, 1]) * dims[2] + keyB[:, 2]
    order = np.argsort(flat)
    flat_s, Bs = flat[order], B[order]
    starts = np.searchsorted(flat_s, np.arange(flat_s[-1] + 2))

    keyA = np.floor((A - lo) / cell).astype(np.int64)
    out = np.full(len(A), np.inf)
    who = np.full(len(A), -1, np.int64)
    for dz in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                k = keyA + (dx, dy, dz)
                ok = np.all((k >= 0) & (k < dims), axis=1)
                if not ok.any(): continue
                f = (k[ok, 0] * dims[1] + k[ok, 1]) * dims[2] + k[ok, 2]
                f = np.clip(f, 0, len(starts) - 2)
                s, e = starts[f], starts[f + 1]
                idx = np.where(ok)[0]
                for i, a, b in zip(idx, s, e):
                    if a == b: continue
                    dd = np.linalg.norm(Bs[a:b] - A[i], axis=1)
                    j = dd.argmin()
                    if dd[j] < out[i]:
                        out[i] = dd[j]
                        who[i] = order[a + j]        # index back into the original B
    return (out, who) if want_index else out

def _selfcheck():
    rng = np.random.default_rng(0)
    P = rng.uniform(0, 100, (4000, 3)); P[:, 2] = 0.0
    Q = P.copy(); Q[:, 2] = 7.0                       # a plane exactly 7 away
    d = nn_dist(P, Q, cell=12.0)
    assert np.allclose(d, 7.0, atol=0.6), (d.min(), d.max())
    Q2 = np.vstack([Q, P + [0, 0, 0.01]])             # add a coincident sheet
    d2 = nn_dist(P, Q2, cell=12.0)
    assert d2.max() < 0.2, d2.max()
    # the index variant must point at the actual nearest row of B
    d3, who = nn_dist(P[:200], Q, cell=12.0, want_index=True)
    assert np.allclose(np.linalg.norm(Q[who] - P[:200], axis=1), d3), "index mismatch"
    print("nn_dist self-check ok (parallel planes -> 7.0, coincident -> ~0)")

if __name__ == "__main__":
    _selfcheck()
