#!/usr/bin/env python3
"""Normalise the three S-LDSC result files (UKB / FinnGen / MVP) into one small
JSON baked into the Space image (gwas_ldsc.json) for the GWAS/Heritability module.

The three files have different (and, for UKB, mis-named) headers, so we index by
COLUMN POSITION. We compute a one-sided coefficient p-value from the z-score and
enrich sample size from the GWAS_annotation files where available.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import csv, json, glob, math, os
from scipy import stats

LDSC = str(config.LDSC_DIR)
ANNO = str(config.GWAS_ANNO_DIR)
OUT  = str(config.GWAS_LDSC_JSON)

# per-source column positions: (file glob, z, coef, se, name, category, id, n_or_None)
SRC = {
    "UKB":     ("*UKB*",     10, 8, 9, 11, 25, 1, None),
    "FinnGen": ("*FinnGen*", 10, 8, 9, 11, 17, 1, None),
    "MVP":     ("*MVP*",     10, 8, 9, 15, 14, 1, 16),
}


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None


# FinnGen sample size from annotation: phenocode -> num_cases + num_controls
finn_n = {}
try:
    with open(ANNO + "/FinnGen_annotation.txt") as fh:
        r = csv.reader(fh, delimiter="\t"); h = next(r)
        for row in r:
            try:
                nc = fnum(row[4]) or 0; ncon = fnum(row[6]) or 0
                finn_n[row[1]] = int(nc + ncon)
            except IndexError:
                pass
except FileNotFoundError:
    pass

out = {"version": "2026-06-29", "metric": "S-LDSC coefficient z-score", "sources": {}}
for src, (g, zi, ci, sei, ni, cati, idi, nci) in SRC.items():
    path = glob.glob(LDSC + "/" + g)[0]
    recs = []
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t"); next(r)
        for row in r:
            if max(zi, ni, cati, idi) >= len(row):
                continue
            z = fnum(row[zi])
            if z is None:
                continue
            n = None
            if nci is not None and nci < len(row):
                n = fnum(row[nci]); n = int(n) if n else None
            if src == "FinnGen":
                n = finn_n.get(row[idi])
            recs.append({
                "id": row[idi], "name": row[ni].strip() or row[idi],
                "category": (row[cati].strip() or "Uncategorized"),
                "coef": fnum(row[ci]), "coef_se": fnum(row[sei]),
                "z": round(z, 4),
                "p": float(stats.norm.sf(z)),          # one-sided (coefficient > 0)
                "n": n,
            })
    # de-duplicate by id (keep the one with the larger |z|)
    best = {}
    for d in recs:
        k = d["id"]
        if k not in best or abs(d["z"]) > abs(best[k]["z"]):
            best[k] = d
    recs = sorted(best.values(), key=lambda d: -d["z"])
    out["sources"][src] = recs
    cats = sorted(set(d["category"] for d in recs))
    print(f"{src:8s}: {len(recs)} traits, {len(cats)} categories, z range [{recs[-1]['z']:.2f}, {recs[0]['z']:.2f}]")

json.dump(out, open(OUT, "w"))
print(f"wrote {OUT}  ({os.path.getsize(OUT)/1e6:.2f} MB)")
