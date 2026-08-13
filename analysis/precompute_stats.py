#!/usr/bin/env python3
"""Precompute genome-wide reference statistics for the Regional-Level Analysis module.

The Panel A/B stratified plots need, per model configuration:
  - the genome-wide mean of |Evo2 delta| (the "genome-wide mean" dashed reference line)
  - its sd / n  (used as the population reference for the one-sample tests)

These are stable (only change when the underlying scores change), so we compute them
once from evo2.db and bake the result into the Space image as genome_stats.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, json, math, time

DB = str(config.DB_PATH)
OUT = str(config.GENOME_STATS_JSON)
DATASET_VERSION = "2026-06-18-gnomad-v4.1-joint"

# (column, friendly label) — the 6 model x strand-strategy delta configurations
CONFIGS = [
    ("Evo2_7B_NoRC_Delta",        "7B · noRC"),
    ("Evo2_7B_AvgRC_Delta",       "7B · avgRC"),
    ("Evo2_7B_WeightedRC_Delta",  "7B · wtRC"),
    ("Evo2_40B_NoRC_Delta",       "40B · noRC"),
    ("Evo2_40B_AvgRC_Delta",      "40B · avgRC"),
    ("Evo2_40B_WeightedRC_Delta", "40B · wtRC"),
]

t0 = time.time()
db = sqlite3.connect(DB)
out = {"dataset_version": DATASET_VERSION, "default_config": "Evo2_40B_NoRC_Delta", "configs": {}}
for col, label in CONFIGS:
    # mean(|x|), mean(x^2) -> sd of |x| in a single scan
    m, m2, n = db.execute(
        f'SELECT avg(abs("{col}")), avg(abs("{col}")*abs("{col}")), count("{col}") FROM variants'
    ).fetchone()
    sd = math.sqrt(max(0.0, (m2 or 0) - (m or 0) ** 2))
    out["configs"][col] = {"label": label, "mean_abs": m, "sd_abs": sd, "n": n}
    print(f"  {label:12s} mean|Δ|={m:.4f}  sd={sd:.4f}  n={n:,}")
db.close()
json.dump(out, open(OUT, "w"), indent=2)
print(f"wrote {OUT}  ({time.time()-t0:.0f}s)")
