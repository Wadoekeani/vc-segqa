"""Figure for the cross-resolution check: xres_results.json -> xres-agreement.png"""
import json, re, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "xres_results.json"))
out = sys.argv[2] if len(sys.argv) > 2 else "xres-agreement.png"
VX = 45.532
CONTROL, FINE24, FINE11 = 5.4, 3.8, 0.34      # xres.py / xres_control.py

def wind(s):
    m = re.search(r"w(\d{3})-(\d{3})", s)
    return (int(m.group(1)) + int(m.group(2))) / 2 if m else None

fig, ax = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw=dict(width_ratios=[1.35, 1]))

w = sorted(((wind(r["seg"]), r) for r in rows if wind(r["seg"])), key=lambda t: t[0])
x = [a for a, _ in w]
ax[0].plot(x, [r["p90"] for _, r in w], "^--", ms=4, lw=.8, color="tab:orange", label="p90")
ax[0].plot(x, [r["p50"] for _, r in w], "o-", ms=4, lw=1.1, color="tab:blue", label="p50")
for y, t, c in [(CONTROL, "measurement noise floor (5.4 µm)", "tab:green"),
                (FINE24, "2.4 vs 7.91 µm, same method (3.8 µm)", "tab:purple")]:
    ax[0].axhline(y, color=c, lw=1, ls=":")
    ax[0].text(x[0], y * 1.12, t, fontsize=7, color=c)
ax[0].axhline(VX, color="grey", lw=.8, ls="--", label="one 45.5 µm voxel")
ax[0].set_yscale("log")
ax[0].set_xlabel("winding number (midpoint of the range in the segment name)")
ax[0].set_ylabel("45.5 µm vs 7.91 µm disagreement (µm)")
ax[0].set_title("the coarse mesh drifts further from the fine one, outward")
ax[0].legend(fontsize=8, loc="lower right")

def grp(s):
    return "w-named" if re.search(r"w\d{3}-\d{3}", s) else ("5753_* family" if "5753_" in s else "everything else")
order = ["w-named", "everything else", "5753_* family"]
data = [[r["p50"] for r in rows if grp(r["seg"]) == g] for g in order]
for i, (g, d) in enumerate(zip(order, data)):
    ax[1].scatter(np.full(len(d), i) + np.random.default_rng(0).uniform(-.12, .12, len(d)),
                  d, s=22, alpha=.8, color=("tab:red" if "5753" in g else "tab:blue"))
    ax[1].plot([i - .28, i + .28], [np.median(d)] * 2, color="k", lw=1.5)
    ax[1].text(i, max(d) + 3, f"n={len(d)}", ha="center", fontsize=8)
ax[1].axhline(VX, color="grey", lw=.8, ls="--")
ax[1].text(2.45, VX + 1.5, "one 45.5 µm voxel", fontsize=7, color="grey", ha="right")
ax[1].set_xticks(range(3)); ax[1].set_xticklabels(order, fontsize=9)
ax[1].set_ylabel("disagreement p50 (µm)"); ax[1].set_ylim(0, 100)
ax[1].set_title("nine segments are twice as bad as the other 46")

fig.tight_layout(); fig.savefig(out, dpi=120)
print("wrote", out)
