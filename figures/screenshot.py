#!/usr/bin/env python3
"""High-res screenshots of the live Datasette pages for the NAR interface figure."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

from playwright.sync_api import sync_playwright
import os

BASE = config.LOCAL_SERVER
OUT = str(config.SCREENS_DIR)
os.makedirs(OUT, exist_ok=True)

# (name, path, full_page)
SHOTS = [
    ("landing",    "/",                                          True),
    ("about",      "/about",                                     True),
    ("gene_apoe",  "/evo2/gene_top_impact?gene=APOE",            False),
    ("faceted",    "/evo2/variants?_facet=Functional_Annotation&_facet=Exonic_Function", False),
    ("variant",    "/evo2/lookup_by_rsid?rsid=rs429358",         False),
    ("search",     "/evo2/variants?_search=APOE",                False),
]

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={"width": 1280, "height": 900}, device_scale_factor=2)
    for name, path, full in SHOTS:
        page.goto(BASE + path, wait_until="networkidle", timeout=30000)
        page.wait_for_timeout(900)
        page.screenshot(path=f"{OUT}/{name}.png", full_page=full)
        print("shot", name, "->", f"{OUT}/{name}.png")
    browser.close()
print("done")
