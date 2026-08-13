#!/usr/bin/env python3
"""Graphical abstract (5:2) for the Evo2 Variant Effect Database."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.colors import LinearSegmentedColormap

D = json.load(open(str(config.CASE_DATA_JSON)))
INK="#0F172A"; ACCENT="#2563EB"; MUTED="#475569"; SLATE="#64748B"; BORDER="#CBD5E1"
EVO = LinearSegmentedColormap.from_list("evo", ["#A32D2D","#E24B4A","#F09595","#C9CCD6","#97C459","#639922"])
plt.rcParams.update({"font.family":"DejaVu Sans"})

fig, ax = plt.subplots(figsize=(12.5, 5.0))   # 5:2
ax.set_xlim(0, 125); ax.set_ylim(0, 50); ax.axis("off")

# ---- DNA double helix (left -> centre), rungs recoloured by diverging scale toward the centre
xs = np.linspace(6, 74, 420)
ph = (xs - 6) / 3.4
y1 = 30 + 6.5*np.sin(ph); y2 = 30 + 6.5*np.sin(ph + np.pi)
ax.plot(xs, y1, color=SLATE, lw=2.6, solid_capstyle="round")
ax.plot(xs, y2, color=ACCENT, lw=2.6, solid_capstyle="round")
for i in range(0, len(xs), 12):
    t = (xs[i] - 6) / (74 - 6)            # 0..1 along helix
    col = SLATE if t < 0.45 else EVO(np.clip((t-0.45)/0.55, 0, 1))
    ax.plot([xs[i], xs[i]], [y1[i], y2[i]], color=col, lw=2.0, alpha=0.9)

# ---- model chips
for j,(lab,yy) in enumerate([("Evo2-7B", 12.5), ("Evo2-40B", 6.0)]):
    ax.add_patch(FancyBboxPatch((7, yy), 18, 4.6, boxstyle="round,pad=0.3,rounding_size=1.2",
                 fc="#EEF2FB", ec=ACCENT, lw=1.3))
    ax.text(16, yy+2.3, lab, ha="center", va="center", fontsize=11, fontweight="bold", color=ACCENT)
ax.text(16, 18.5, "DNA foundation model", ha="center", fontsize=8.5, color=MUTED)

# ---- centre: equation + mini delta histogram
ax.text(54, 44, "Δ = log P(ALT) − log P(REF)", fontsize=12.5, fontweight="bold", color=INK, ha="center")
d40 = np.array(D["hist"]["d40"]); d40 = d40[np.abs(d40) < 120]
counts, edges = np.histogram(d40, bins=30)
hx = (edges[:-1]+edges[1:])/2; hw = edges[1]-edges[0]
hxn = 42 + (hx - hx.min())/(hx.max()-hx.min())*30      # map to x 42..72
hh = counts/counts.max()*12                            # height up to 12
for cx, h, raw in zip(hxn, hh, hx):
    ax.add_patch(Rectangle((cx-0.45, 6), 0.9, h, color=EVO(np.clip((raw+60)/120,0,1)), ec="none"))
ax.plot([42,72],[6,6], color=BORDER, lw=1.0)
ax.text(57, 2.6, "heavy-tailed Δ distribution", ha="center", fontsize=8.5, color=MUTED)

# ---- arrow to outputs
ax.annotate("", xy=(83, 30), xytext=(76, 30),
            arrowprops=dict(arrowstyle="-|>", lw=2.2, color=SLATE))

# ---- right: access icons + headline
bx = 88
# browser
ax.add_patch(FancyBboxPatch((bx, 36), 9, 6, boxstyle="round,pad=0.2,rounding_size=0.8", fc="white", ec=ACCENT, lw=1.3))
ax.add_patch(Rectangle((bx, 40), 9, 2, color=ACCENT)); ax.text(bx+4.5, 37.6, "browse", ha="center", fontsize=7.5, color=INK)
# SQL
ax.add_patch(FancyBboxPatch((bx+11, 36), 9, 6, boxstyle="round,pad=0.2,rounding_size=0.8", fc="white", ec=ACCENT, lw=1.3))
ax.text(bx+15.5, 39, "SQL", ha="center", va="center", fontsize=9, fontweight="bold", color=ACCENT)
# API
ax.add_patch(FancyBboxPatch((bx+22, 36), 9, 6, boxstyle="round,pad=0.2,rounding_size=0.8", fc="white", ec=ACCENT, lw=1.3))
ax.text(bx+26.5, 39, "{ } API", ha="center", va="center", fontsize=8, fontweight="bold", color=ACCENT)

ax.text(104, 27, "6,475,578", ha="center", fontsize=20, fontweight="bold", color=INK)
ax.text(104, 22, "common human variants", ha="center", fontsize=10.5, color=INK)
ax.text(104, 17.5, "scored zero-shot by a", ha="center", fontsize=10.5, color=MUTED)
ax.text(104, 13.5, "DNA foundation model", ha="center", fontsize=10.5, color=MUTED)

# ---- title strip
ax.text(2, 48.2, "Evo2 Variant Effect Database", fontsize=14, fontweight="bold", color=INK)

plt.tight_layout()
for ext, kw in [("pdf", {}), ("png", {"dpi": 300})]:
    fig.savefig(f"{config.FIGURE_DIR}/Graphical_abstract.{ext}", bbox_inches="tight", **kw)
print("graphical abstract saved")
