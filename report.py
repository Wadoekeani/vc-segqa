"""Summarise a corpus scan: worst offenders, per-scroll rates, winding trend."""
import re, sys, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = json.load(open(sys.argv[1]))
out = sys.argv[2] if len(sys.argv) > 2 else "scan-summary.png"
scrolls = sorted({r["scroll"] for r in rows})

print(f"{'scroll':14s} {'segs':>5} {'verts':>11} {'outside':>9} {'>5%':>5} {'worst':>8}")
per = {}
for s in scrolls:
    rs = [r for r in rows if r["scroll"] == s]
    tot = sum(r["verts"] for r in rs)
    bad = sum(r["verts"] * r["outside"] for r in rs)
    per[s] = (len(rs), tot, bad/max(tot, 1), sum(r["outside"] > .05 for r in rs),
              max(r["outside"] for r in rs))
    print(f"{s:14s} {len(rs):>5} {tot:>11,} {per[s][2]:>8.2%} {per[s][3]:>5} {per[s][4]:>7.1%}")
tot = sum(r["verts"] for r in rows)
bad = sum(r["verts"]*r["outside"] for r in rows)
print(f"\n{'ALL':14s} {len(rows):>5} {tot:>11,} {bad/tot:>8.2%} "
      f"{sum(r['outside'] > .05 for r in rows):>5}")

print(f"\nworst 15 segments")
for r in sorted(rows, key=lambda r: -r["outside"])[:15]:
    print(f"  {r['scroll']:13s} {r['segment'][:46]:46s} {r['verts']:>8,} {r['outside']:>7.2%}")

# winding trend: segment names carry wNNN
fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
ax[0].barh(scrolls, [per[s][2]*100 for s in scrolls], color="tab:red", alpha=.8)
ax[0].set_xlabel("% of vertices outside the scroll mask"); ax[0].set_title("per scroll")
for s in scrolls:
    pts = []
    for r in rows:
        if r["scroll"] != s: continue
        m = re.search(r"w(\d{3})", r["segment"])
        if m: pts.append((int(m.group(1)), r["outside"]*100))
    if len(pts) > 3:
        pts.sort(); a = np.array(pts)
        ax[1].plot(a[:, 0], a[:, 1], ".-", ms=4, lw=.9, label=s, alpha=.85)
ax[1].set_xlabel("winding number (from the segment name)")
ax[1].set_ylabel("% outside mask"); ax[1].set_title("outer windings drift off the object")
ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(out, dpi=120); print("\nwrote", out)
