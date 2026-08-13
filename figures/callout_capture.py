#!/usr/bin/env python3
"""Capture cropped panel screenshots + precise DOM element boxes for Fig2 callouts."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

from playwright.sync_api import sync_playwright
import json, os

BASE=config.LOCAL_SERVER; OUT=str(config.SCREENS_DIR)
SCALE=2; os.makedirs(OUT, exist_ok=True)

# panel -> (path, [(css_selector, label)])
PANELS = [
    ("p_landing", "/", [(".evo2-search","Four entry points: rsID · gene · region · free-text"),
                        (".evo2-legend-bar","Diverging Δ scale: red disfavored → green tolerated")]),
    ("p_viewer", "/viewer?chr=19&start=43900000&end=45900000", [("#landscape","Interactive Δ landscape (APOE locus): position vs Evo2-40B Δ, red→green colour scale")]),
    ("p_faceted", "/evo2/variants?_facet=Functional_Annotation", [(".facet-results","Facet by functional region, with live counts")]),
    ("p_about", "/about", [(".evo2-doc table","2 models × 3 strand strategies → 18 model-score columns")]),
]

result=[]
with sync_playwright() as p:
    b=p.chromium.launch(); pg=b.new_page(viewport={"width":1180,"height":900}, device_scale_factor=SCALE)
    for name, path, targets in PANELS:
        pg.goto(BASE+path, wait_until="networkidle", timeout=30000); pg.wait_for_timeout(700)
        if name=="p_viewer":
            try: pg.wait_for_selector("#landscape canvas, #landscape svg", timeout=15000)
            except Exception: pass
            pg.wait_for_timeout(1800)
        full=f"{OUT}/{name}_full.png"; pg.screenshot(path=full, full_page=True)
        boxes=[]
        for sel,label in targets:
            loc=pg.locator(sel).first
            if loc.count():
                bb=loc.bounding_box()
                if bb: boxes.append({"css":bb,"label":label})
        if not boxes:
            result.append({"name":name,"img":full,"crop":[0,0],"boxes":[]}); continue
        # crop window (css px) framing all target boxes
        top=max(0, min(b["css"]["y"] for b in boxes)-45)
        bot=max(b["css"]["y"]+min(b["css"]["height"],480) for b in boxes)+55
        # image px
        import PIL.Image as I
        im=I.open(full); W,H=im.size
        ct,cb=int(top*SCALE), min(H,int(bot*SCALE))
        crop=im.crop((0,ct,W,cb)); cpath=f"{OUT}/{name}_crop.png"; crop.save(cpath)
        pix=[{"x":b["css"]["x"]*SCALE,"y":b["css"]["y"]*SCALE-ct,
              "w":b["css"]["width"]*SCALE,"h":b["css"]["height"]*SCALE,"label":b["label"]} for b in boxes]
        result.append({"name":name,"img":cpath,"boxes":pix})
        print(name,"boxes:",len(pix))
    b.close()
json.dump(result, open(f"{OUT}/fig2_panels.json","w"))
print("captured")
