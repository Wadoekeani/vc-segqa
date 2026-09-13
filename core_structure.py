"""What does the 2.4 um image say about Scroll 1's untraced core?

The innermost traced sheet is w010; nothing exists inside it. Two ways to read
that - nobody has tried, or the image there cannot support a trace. This asks
the image directly, at twelve heights, with the traced band just outside w010 as
the control on the same slice:

  coverage   fraction of the area whose gradient energy clears a global
             threshold - how much of it has any structure at all
  coherence  |mean unit double-angle gradient vector| over a 96 um window -
             do those structures line up into sheets, or point every which way

Coherence deliberately uses no centre. A concentricity measure was tried first
and swung from 0.41 to 0.87 when the assumed axis moved 60 voxels, which at
small radius it always might. The two measures kept here only use the axis to
draw the bands; __main__ ends with two sensitivity blocks - a different
"structure" threshold, and the axis shifted by 60 voxels - so the size of that
dependence is printed rather than asserted.

The traced band is w010's inner radius out to 1.8x that radius, capped at the
window edge. Where w010 already sits at the window edge (two of the twelve
heights) there is no band and the control is NaN; those heights count in the
core numbers only.

Needs the 45.5 um w010-027 meshes (for the axis) and the 45.5 um and 2.4 um
transform.json files (xres_fetch.py). Pulls ~2 MB per height from the 2.4 um
volume over HTTP Range; nothing else.
"""
import json, glob, os, sys, numpy as np
import xres
from coverage import load_seg
from axis import fit_axis, radius
from window import plane_window

VOL = "20260411134726-2.400um-0.2m-78keV-masked.zarr"
HALF = 620                    # window half-size in 2.4 um voxels (1.5 mm)
BIN = 2                       # 2x2 block mean before measuring: 4.8 um/px
WINDOW_UM = 96.0              # coherence window; sheets here are 180-400 um apart
GATE_PCT = 60                 # "has structure" = gradient energy above this percentile
CACHE = xres.CACHE


def _blur(a, sig):
    r = int(3 * sig)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sig) ** 2); k /= k.sum()
    out = np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 0, a)
    return np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 1, out)


def core_geometry(heights):
    """For each 45.5 um z: the scroll axis in 2.4 um voxels and w010's inner
    radius there. Axis from the w010-027 meshes, mapped 45.5 -> fixed -> 2.4."""
    P = np.concatenate([load_seg(d) for d in sorted(glob.glob(f"{CACHE}/PHercParis4/*w010-027*/"))])
    zc, cen = fit_axis(P); r = radius(P, zc, cen)
    M = np.array(json.load(open(f"{xres.TDIR}/{xres.VOLS['2.4']}.json"))["transformation_matrix"])
    inv = np.linalg.inv(M[:, :3])
    to24 = lambda F: (F - M[:, 3]) @ inv.T
    f45 = xres.to_fixed(xres.VOLS["45.5"])
    out = []
    for z45 in heights:
        near = np.abs(P[:, 2] - z45) < 60
        if near.sum() < 300: continue
        k = np.argmin(abs(zc - z45))
        rin = np.percentile(r[near], 1) * 45.532 / 2.4
        c = to24(f45(np.array([[cen[k, 0], cen[k, 1], float(z45)]])))[0]
        out.append(dict(z45=z45, cx=int(round(c[0])), cy=int(round(c[1])), z=int(round(c[2])), rin=rin))
    return out


def fetch(g, cache_dir):
    f = f"{cache_dir}/core-z{g['z45']}.npy"
    if os.path.exists(f):
        return np.load(f), tuple(np.load(f.replace(".npy", "-org.npy")))
    img, y0, x0, _ = plane_window("PHercParis4", VOL, "0", g["z"], g["cy"], g["cx"], HALF)
    np.save(f, img); np.save(f.replace(".npy", "-org.npy"), np.array([y0, x0]))
    return img, (y0, x0)


def measure(img, y0, x0, cx, cy, rin, gate_pct=GATE_PCT):
    a = img.astype(np.float64)
    h, w = (a.shape[0] // BIN) * BIN, (a.shape[1] // BIN) * BIN
    a = a[:h, :w].reshape(h // BIN, BIN, w // BIN, BIN).mean((1, 3))
    b = _blur(a, 1.0); gy, gx = np.gradient(b)
    mag = np.hypot(gx, gy); th = np.arctan2(gy, gx)
    gate = mag > np.percentile(mag, gate_pct)
    wgt = gate.astype(float); sig = WINDOW_UM / (2.4 * BIN)
    coh = np.hypot(_blur(np.cos(2 * th) * wgt, sig), _blur(np.sin(2 * th) * wgt, sig)) \
        / np.maximum(_blur(wgt, sig), 1e-9)
    yy, xx = np.mgrid[0:a.shape[0], 0:a.shape[1]]
    R = np.hypot(x0 + xx * BIN - cx, y0 + yy * BIN - cy)

    def band(r0, r1):
        m = (R >= r0) & (R < r1); g = m & gate
        if m.sum() < 500 or g.sum() < 200: return np.nan, np.nan
        return g.sum() / m.sum(), float(np.median(coh[g]))
    return band(30, rin), band(rin, min(rin * 1.8, HALF))     # untraced core | traced band


if __name__ == "__main__":
    cache_dir = f"{CACHE}/_core"; os.makedirs(cache_dir, exist_ok=True)
    rows = []
    for g in core_geometry(range(800, 3900, 280)):
        img, (y0, x0) = fetch(g, cache_dir)
        (cov_i, coh_i), (cov_o, coh_o) = measure(img, y0, x0, g["cx"], g["cy"], g["rin"])
        rows.append(dict(z45=g["z45"], z24=g["z"], rin_um=round(g["rin"] * 2.4),
                         cov_core=cov_i, coh_core=coh_i, cov_traced=cov_o, coh_traced=coh_o))
        print(f"z45={g['z45']:4} r_in={g['rin']*2.4/1000:.2f} mm | core cov={cov_i:6.1%} coh={coh_i:.3f} "
              f"| traced cov={cov_o:6.1%} coh={coh_o:.3f}", flush=True)
    json.dump(rows, open("core_structure.json", "w"), indent=1)

    z = np.array([r["z45"] for r in rows])
    cc, ct = np.array([r["cov_core"] for r in rows]), np.array([r["cov_traced"] for r in rows])
    hc, ht = np.array([r["coh_core"] for r in rows]), np.array([r["coh_traced"] for r in rows])
    ok = np.isfinite(ct)
    dcoh = (hc - ht)[ok]
    print(f"\ncoherence, core minus traced band: median {np.median(dcoh):+.3f}  "
          f"range {dcoh.min():+.3f}..{dcoh.max():+.3f}  (n={ok.sum()})")
    lo = ok & (z <= 2200)
    print(f"coverage, lower half (z45<=2200): core {np.mean(cc[lo]):.1%} vs traced {np.mean(ct[lo]):.1%}")
    print(f"coverage, upper half (z45> 2200): core {np.mean(cc[ok & ~lo]):.1%} vs traced {np.mean(ct[ok & ~lo]):.1%}")

    # sensitivity 1: the conclusion has to survive a different "this is papyrus"
    # threshold, or it is a threshold artefact - for coverage as well as coherence
    geoms = core_geometry([r["z45"] for r in rows])
    for pct in (40, 80):
        dcoh, dcov_lo, dcov_hi = [], [], []
        for g in geoms:
            img, (y0, x0) = fetch(g, cache_dir)
            (ci, a), (co, b) = measure(img, y0, x0, g["cx"], g["cy"], g["rin"], pct)
            if not (np.isfinite(a) and np.isfinite(b)): continue
            dcoh.append(a - b); (dcov_lo if g["z45"] <= 2200 else dcov_hi).append(ci - co)
        print(f"  gate at {pct}th pct: coherence core-traced median {np.median(dcoh):+.3f} | "
              f"coverage core-traced lower {100*np.median(dcov_lo):+.1f} pts, upper {100*np.median(dcov_hi):+.1f} pts")
    # sensitivity 2: the bands come from a fitted axis; shift it 60 voxels (144 um)
    # in four directions and see how far the paired differences move
    dcoh, dcov = [], []
    for g in geoms:
        img, (y0, x0) = fetch(g, cache_dir)
        (ci0, a0), (co0, b0) = measure(img, y0, x0, g["cx"], g["cy"], g["rin"])
        if not (np.isfinite(b0)): continue
        for dx, dy in ((60, 0), (-60, 0), (0, 60), (0, -60)):
            (ci, a), (co, b) = measure(img, y0, x0, g["cx"] + dx, g["cy"] + dy, g["rin"])
            if np.isfinite(b): dcoh.append((a - b) - (a0 - b0)); dcov.append((ci - co) - (ci0 - co0))
    print(f"  axis shifted 60 vx: paired coherence difference moves by median {np.median(np.abs(dcoh)):.3f}, "
          f"coverage difference by median {100*np.median(np.abs(dcov)):.1f} pts  (n={len(dcoh)})")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    mm = z * 45.532 / 1000
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    ax[0].plot(mm, cc * 100, "o-", color="tab:red", label="untraced core")
    ax[0].plot(mm, ct * 100, "s--", color="tab:grey", label="traced band (w010+)")
    ax[0].set_ylabel("% of area with structure"); ax[0].set_title("how much papyrus")
    ax[1].plot(mm, hc, "o-", color="tab:red", label="untraced core")
    ax[1].plot(mm, ht, "s--", color="tab:grey", label="traced band (w010+)")
    ax[1].set_ylabel("orientation coherence"); ax[1].set_title("does it line up into sheets")
    ax[2].plot(mm, [r["rin_um"] / 1000 for r in rows], "o-", color="tab:blue")
    ax[2].set_ylabel("radius of the innermost traced sheet (mm)"); ax[2].set_title("how much is left untraced")
    for a in ax: a.set_xlabel("height along the scroll axis (mm)"); a.grid(alpha=.25)
    ax[0].legend(fontsize=8); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig("core-structure.png", dpi=110)
    print("wrote core_structure.json, core-structure.png")
