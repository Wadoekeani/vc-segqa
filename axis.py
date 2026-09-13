"""Scroll axis as a function of height, then radius measured against it.

The scroll is bent: a single global centre puts points near the axis at a
meaningless radius. Estimate the umbilicus per z-slab from the centroid of
every traced vertex in that slab - the traced body spans many windings and the
full angular range, so its centroid tracks the core.
"""
import numpy as np

def fit_axis(pts, step=50.0):
    z = pts[:, 2]
    z0, z1 = z.min(), z.max()
    edges = np.arange(z0, z1 + step, step)
    idx = np.clip(np.searchsorted(edges, z) - 1, 0, len(edges) - 2)
    cen = np.full((len(edges) - 1, 2), np.nan)
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() >= 200: cen[b] = pts[m, :2].mean(0)
    mid = (edges[:-1] + edges[1:]) / 2
    ok = np.isfinite(cen).all(1)
    return mid[ok], cen[ok]                       # (nb,), (nb, 2)

def radius(pts, zc, cen):
    cx = np.interp(pts[:, 2], zc, cen[:, 0])
    cy = np.interp(pts[:, 2], zc, cen[:, 1])
    return np.hypot(pts[:, 0] - cx, pts[:, 1] - cy)
