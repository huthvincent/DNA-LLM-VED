#!/usr/bin/env python3
"""Figure 2 — web interface, 2x2 cropped screenshots with precise DOM-anchored callouts."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch, Circle

S=str(config.SCREENS_DIR)
P=json.load(open(f"{S}/fig2_panels.json"))
INK="#0F172A"; CALL="#C0392B"
titles={"p_landing":("a","Landing — search entry points & Δ legend"),
        "p_viewer":("b","Region viewer — interactive Δ landscape (APOE locus)"),
        "p_faceted":("c","Faceted browse of all 6.48M variants"),
        "p_about":("d","About / feature scheme")}
order=["p_landing","p_viewer","p_faceted","p_about"]
by={p["name"]:p for p in P}

fig,axs=plt.subplots(2,2,figsize=(13,10.2)); plt.rcParams.update({"font.family":"DejaVu Sans"})
num=0; legend=[]
for ax,key in zip(axs.flat,order):
    panel=by[key]; img=mpimg.imread(panel["img"]); H,W=img.shape[0],img.shape[1]
    ax.imshow(img,extent=[0,W,H,0]); ax.set_xlim(0,W); ax.set_ylim(H,0)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_edgecolor("#CBD5E1"); sp.set_linewidth(1.0)
    lett,title=titles[key]
    ax.set_title(title,fontsize=10,color=INK,loc="left",pad=6)
    ax.text(-0.01,1.11,lett,transform=ax.transAxes,fontsize=15,fontweight="bold",color=INK,ha="right",va="bottom")
    for box in panel["boxes"]:
        num+=1
        bh=min(box["h"], H-box["y"]-6)
        ax.add_patch(FancyBboxPatch((box["x"],box["y"]),box["w"],bh,
                     boxstyle="round,pad=2,rounding_size=6",fill=False,edgecolor=CALL,linewidth=2.2))
        cx,cy=box["x"],box["y"]
        ax.add_patch(Circle((cx,cy),16,color=CALL,zorder=5))
        ax.text(cx,cy,str(num),color="white",fontsize=10,fontweight="bold",ha="center",va="center",zorder=6)
        legend.append(f"{num}  {box['label']}")

fig.subplots_adjust(bottom=0.13,hspace=0.16,wspace=0.08)
fig.text(0.5,0.065," · ".join([]) ,ha="center")  # spacer
for i,l in enumerate(legend):
    fig.text(0.07, 0.085 - i*0.018, l, fontsize=9, color=INK, ha="left", va="top",
             fontfamily="DejaVu Sans")
for ext,kw in [("pdf",{}),("png",{"dpi":200})]:
    fig.savefig(f"{config.FIGURE_DIR}/Fig2_interface.{ext}",bbox_inches="tight",**kw)
print("Fig2 saved with",num,"callouts")
