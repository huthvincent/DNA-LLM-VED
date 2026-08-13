#!/usr/bin/env python3
"""Fig2 — each panel as a standalone high-res (3x) screenshot."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

from playwright.sync_api import sync_playwright

BASE=config.LOCAL_SERVER; OUT=str(config.FIGURE_DIR / "Fig2"); W=1180

def clip(pg,name,x,y,w,h): pg.screenshot(path=f"{OUT}/{name}.png",clip={"x":x,"y":max(0,y),"width":w,"height":h})

with sync_playwright() as p:
    b=p.chromium.launch()
    pg=b.new_page(viewport={"width":W,"height":2300},device_scale_factor=3)

    # a — landing entry points + Δ legend  (top -> bottom of the "what am I looking at" note)
    pg.goto(BASE+"/",wait_until="networkidle"); pg.wait_for_timeout(600)
    nb=pg.locator(".evo2-section").first.bounding_box()  # "What am I looking at?" + Δ legend
    clip(pg,"a_landing_entrypoints",0,0,W,nb["y"]+nb["height"]+24)
    print("a done")

    # b — region viewer (form + the 3 interactive charts)
    pg.goto(BASE+"/viewer?chr=19&start=43900000&end=45900000",wait_until="networkidle",timeout=40000)
    try: pg.wait_for_selector("#landscape canvas, #landscape svg",timeout=20000)
    except Exception: pass
    pg.wait_for_timeout(2500)
    f=pg.locator(".evo2-viewer-form").bounding_box(); h=pg.locator("#hist").bounding_box()
    clip(pg,"b_region_viewer",0,f["y"]-16,W,(h["y"]+h["height"]+24)-(f["y"]-16))
    print("b done")

    # c — faceted browse (facet sidebar + table top)
    pg.goto(BASE+"/evo2/variants?_facet=Functional_Annotation&_facet=Exonic_Function",wait_until="networkidle"); pg.wait_for_timeout(700)
    clip(pg,"c_faceted_browse",0,0,W,1140)
    print("c done")

    # d — about / feature scheme (the models x strategies + 18-column tables)
    pg.goto(BASE+"/about",wait_until="networkidle"); pg.wait_for_timeout(600)
    hm=pg.locator("h2:has-text('Two models')").bounding_box()
    clip(pg,"d_about_feature_scheme",0,hm["y"]-14,W,660)
    print("d done")

    b.close()
print("Fig2 panels done")
