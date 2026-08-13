#!/usr/bin/env python3
"""Figure 3 -- APOE-locus case study (4 panels, all drawn from ONE region).

The whole figure characterises a single locus: the APOE region +/- 1 Mb on chr19
(~4.8k common variants). Panels a-d all describe that same set of variants:
(a) Delta by functional region, (b) Delta by coding consequence (only the classes
that occur in common variants), (c) score concordance within the locus, and
(d) the variant landscape with the APOE 167-site core highlighted."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

D = json.load(open(str(config.CASE_DATA_JSON)))
ACCENT="#2563EB"; INK="#0F172A"; MUTED="#475569"; GRID="#E2E8F0"
EVO = LinearSegmentedColormap.from_list("evo_delta", ["#A32D2D","#E24B4A","#F09595","#C9CCD6","#97C459","#639922"])
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.edgecolor":"#94A3B8",
                     "axes.linewidth":0.8,"axes.titlesize":11,"axes.titleweight":"bold"})

RG = D["region"]; LOCUS = f"chr{RG['chr']}:{RG['start']/1e6:.1f}-{RG['end']/1e6:.1f} Mb"


def letter(ax,s): ax.text(-0.14,1.05,s,transform=ax.transAxes,fontsize=14,fontweight="bold",color=INK,va="top")


def boxes(ax, groups, labels, title, ylab="Evo2-40B Δ score"):
    data=[np.array(g) for g in groups]; means=[float(np.mean(g)) if len(g) else 0 for g in groups]
    vmax=max(1.0, max(abs(m) for m in means)); norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    bp=ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.62,
                  medianprops={"color":INK,"linewidth":1.3})
    for patch,m in zip(bp["boxes"],means):
        patch.set(facecolor=EVO(norm(m)), edgecolor="#334155", linewidth=0.8)
    ax.axhline(0,color=MUTED,lw=0.8,ls="--"); ax.set_yscale("symlog",linthresh=10)
    ax.set_xticks(range(1,len(labels)+1)); ax.set_xticklabels(labels,rotation=25,ha="right",fontsize=8.5)
    ax.set_ylabel(ylab); ax.set_title(title); ax.grid(axis="y",color=GRID,lw=0.6); ax.set_axisbelow(True)
    return means


fig,axs=plt.subplots(2,2,figsize=(12,9.4))

# A. delta by functional region (within the locus)
ax=axs[0,0]
forder=[f for f in ["exonic","UTR5","UTR3","ncRNA_exonic","upstream","downstream","intronic","ncRNA_intronic","intergenic"] if f in D["delta_by_func"]]
flab=[f.replace("ncRNA_","ncRNA ") for f in forder]
boxes(ax,[D["delta_by_func"][f] for f in forder],flab,"Δ by functional region (APOE locus)")
for i,f in enumerate(forder):
    ax.text(i+1, ax.get_ylim()[1], f"n={len(D['delta_by_func'][f])}", ha="center",va="bottom",fontsize=6.5,color=MUTED)
letter(ax,"a")

# B. delta by coding consequence (within the locus; only classes present in common variants)
ax=axs[0,1]
eorder=[e for e in ["synonymous SNV","nonsynonymous SNV","startloss","stoploss","stopgain"] if e in D["delta_by_exonic"]]
elab={"synonymous SNV":"synonymous","nonsynonymous SNV":"nonsynonymous"}
boxes(ax,[D["delta_by_exonic"][e] for e in eorder],[elab.get(e,e) for e in eorder],"Δ by coding consequence (APOE locus)")
for i,e in enumerate(eorder):
    ax.text(i+1, ax.get_ylim()[1], f"n={len(D['delta_by_exonic'][e])}", ha="center",va="bottom",fontsize=6.5,color=MUTED)
ax.text(0.5,-0.32,"Common variants: only synonymous & nonsynonymous occur\n(protein-truncating changes can't reach 5% frequency)",
        transform=ax.transAxes,ha="center",va="top",fontsize=7,color=MUTED)
letter(ax,"b")

# C. Spearman concordance among all 6 delta scores, computed WITHIN the locus
ax=axs[1,0]
dcols=["Evo2_7B_NoRC_Delta","Evo2_7B_AvgRC_Delta","Evo2_7B_WeightedRC_Delta",
       "Evo2_40B_NoRC_Delta","Evo2_40B_AvgRC_Delta","Evo2_40B_WeightedRC_Delta"]
A=np.array([[r[c] for c in dcols] for r in D["region_scores"]],dtype=float)
R=np.apply_along_axis(lambda c: np.argsort(np.argsort(c)),0,A)
M=np.corrcoef(R,rowvar=False)
slab=["7B noRC","7B avgRC","7B wtRC","40B noRC","40B avgRC","40B wtRC"]
im=ax.imshow(M,cmap="Blues",vmin=0.2,vmax=1.0)
ax.set_xticks(range(6)); ax.set_xticklabels(slab,rotation=45,ha="right",fontsize=8)
ax.set_yticks(range(6)); ax.set_yticklabels(slab,fontsize=8)
for i in range(6):
    for j in range(6):
        ax.text(j,i,f"{M[i,j]:.2f}",ha="center",va="center",fontsize=7.5,
                color="white" if M[i,j]>0.7 else INK)
cb=plt.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.set_label("Spearman ρ",fontsize=8)
ax.set_title(f"Score concordance within the locus (n={len(A)})")
letter(ax,"c")

# D. variant landscape across the locus, with the APOE 167-site core highlighted
ax=axs[1,1]
gc=D["gene_case"]; b=gc["rows"]; chrom=b[0]["CHR"]; core=gc["core"]
bx=np.array([r_["BP"] for r_ in b])/1e6; bd=np.array([r_["d40"] for r_ in b])
vlim=max(abs(np.percentile(bd,2)),abs(np.percentile(bd,98)),1); norm=TwoSlopeNorm(vcenter=0,vmin=-vlim,vmax=vlim)
ax.axvspan(core["start"]/1e6, core["end"]/1e6, color="#2563EB", alpha=0.10, zorder=0)
ax.text((core["start"]+core["end"])/2e6, 0.98, "APOE core", transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=7.5, color=ACCENT)
sc=ax.scatter(bx,bd,c=bd,cmap=EVO,norm=norm,s=18,edgecolors="#334155",linewidths=0.15,alpha=0.85,zorder=3)
ax.axhline(0,color=MUTED,lw=0.8,ls="--"); ax.set_yscale("symlog",linthresh=10)
ax.set_xlabel(f"chr{chrom} position (Mb)"); ax.set_ylabel("Evo2-40B Δ score (symlog)")
ax.set_title(f"Variant landscape across the APOE locus (n={len(b)} common)")
cb=plt.colorbar(sc,ax=ax,fraction=0.046,pad=0.03); cb.set_label("Δ  (disfavored ↔ tolerated)",fontsize=8)
ax.grid(color=GRID,lw=0.6); ax.set_axisbelow(True); letter(ax,"d")

plt.tight_layout(w_pad=2.6,h_pad=3.2)
for ext,kw in [("pdf",{}),("png",{"dpi":300})]:
    fig.savefig(f"{config.FIGURE_DIR}/Fig3_casestudy.{ext}",bbox_inches="tight",**kw)
print(f"Fig3 saved. region={LOCUS} n={len(b)} | cross-scale 7B-40B noRC rho={M[0,3]:.2f}; within-7B noRC-avgRC rho={M[0,1]:.2f}")
