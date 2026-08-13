#!/usr/bin/env python3
"""Figure 1 — Evo2 Variant Effect Database data-processing & access pipeline."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})

INK = "#0F172A"; ACCENT = "#1D4ED8"; BODY = "#1E293B"; MUTED = "#475569"
FILL = "#F1F5F9"; EDGE = "#CBD5E1"; ARROW = "#64748B"; TOPBAR = "#2563EB"

stages = [
    ("1  Input variants", ["gnomAD NFE common variants", "MAF ≥ 5%, call-rate QC",
                            "GRCh38, chr1–22", "6,475,578 sites"]),
    ("2  Evo2 scoring", ["Evo2-7B  &  Evo2-40B", "3 strand strategies:", "noRC / avgRC / weightRC",
                         "ALT & REF log-likelihood → Δ", "18 score columns"]),
    ("3  Annotation", ["ANNOVAR / Ensembl genes", "ENCODE cCRE category",
                       "phastCons / phyloP 100-way", "RepeatMasker",
                       "gnomAD v4.1 joint AF (9 pop.)"]),
    ("4  Database build", ["SQLite · 6.48M rows", "FTS5 full-text + indexes",
                           "Variant_ID = Chr:Pos:Ref:Alt"]),
    ("5  Web resource", ["Datasette app", "search: rsID / gene / region",
                         "faceted browse · SQL · JSON API", "free, no login"]),
]

fig, ax = plt.subplots(figsize=(13.2, 4.3))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

n = len(stages); gap = 2.2
bw = (100 - gap * (n - 1)) / n
y0, bh = 20, 60
for i, (title, items) in enumerate(stages):
    x = i * (bw + gap)
    box = FancyBboxPatch((x, y0), bw, bh, boxstyle="round,pad=0.4,rounding_size=2.2",
                         linewidth=1.1, edgecolor=EDGE, facecolor=FILL, mutation_aspect=0.5)
    ax.add_patch(box)
    ax.add_patch(plt.Rectangle((x + 1.2, y0 + bh - 7), bw - 2.4, 2.0, color=TOPBAR, zorder=3))
    ax.text(x + bw / 2, y0 + bh - 13, title, ha="center", va="center",
            fontsize=12.5, fontweight="bold", color=ACCENT)
    for j, it in enumerate(items):
        ax.text(x + bw / 2, y0 + bh - 21 - j * 7.4, it, ha="center", va="center",
                fontsize=9.3, color=BODY)
    if i < n - 1:
        ax.add_patch(FancyArrowPatch((x + bw + 0.1, y0 + bh / 2), (x + bw + gap - 0.1, y0 + bh / 2),
                     arrowstyle="-|>", mutation_scale=18, linewidth=2.0, color=ARROW))

ax.text(0, 92, "Evo2 Variant Effect Database", fontsize=15, fontweight="bold", color=INK)
ax.text(0, 85.5, "From genome-scale model inference to a free, queryable web resource",
        fontsize=10.5, color=MUTED)

plt.tight_layout()
for ext, kw in [("pdf", {}), ("png", {"dpi": 300})]:
    fig.savefig(f"{config.FIGURE_DIR}/Fig1_architecture.{ext}",
                bbox_inches="tight", **kw)
print("Fig1 saved (pdf + png)")
