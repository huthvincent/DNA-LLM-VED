#!/usr/bin/env python3
"""Add a MAF_tier column to evo2.db by quintile-binning the minor-allele frequency.

Per the July-11 instruction, the current (gnomAD v4.1 joint) file has no MAF_tier,
so we derive it from the Non-Finnish-European allele frequency (the cohort the
variant set was selected on): MAF = min(NFE_AF, 1 - NFE_AF), then Q1..Q5 = quintiles.
Category labels match the reference R code (Q1_0-20% ... Q5_80-100%).

Also writes a small annotated file (MAF_tier_annotation.tsv.gz) for inspection.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, numpy as np, gzip, csv, time

DB = str(config.DB_PATH)
OUT = str(config.MAF_TIER_TSV)
LABELS = ["Q1_0-20%", "Q2_20-40%", "Q3_40-60%", "Q4_60-80%", "Q5_80-100%"]

t0 = time.time()
db = sqlite3.connect(DB)

# 1) minor-allele frequency from NFE AF, and quintile boundaries
maf = np.array([min(a, 1 - a) for (a,) in
                db.execute("SELECT Non_Finnish_European_AF FROM variants WHERE Non_Finnish_European_AF IS NOT NULL")])
q = np.quantile(maf, [0.2, 0.4, 0.6, 0.8])
print(f"MAF quintile cut-points (NFE): {[round(x,5) for x in q]}  (n={len(maf):,}, {time.time()-t0:.0f}s)")

# 2) add column + assign tier via a single CASE update (MAF = min(af, 1-af))
cols = [r[1] for r in db.execute("PRAGMA table_info(variants)")]
if "MAF_tier" not in cols:
    db.execute("ALTER TABLE variants ADD COLUMN MAF_tier TEXT")
mafexpr = "min(Non_Finnish_European_AF, 1-Non_Finnish_European_AF)"
db.execute(f"""
UPDATE variants SET MAF_tier = CASE
  WHEN Non_Finnish_European_AF IS NULL THEN NULL
  WHEN {mafexpr} < {q[0]} THEN '{LABELS[0]}'
  WHEN {mafexpr} < {q[1]} THEN '{LABELS[1]}'
  WHEN {mafexpr} < {q[2]} THEN '{LABELS[2]}'
  WHEN {mafexpr} < {q[3]} THEN '{LABELS[3]}'
  ELSE '{LABELS[4]}' END
""")
db.execute("CREATE INDEX IF NOT EXISTS idx_maf_tier ON variants(MAF_tier)")
db.commit()

print("tier counts:")
for tier, n in db.execute("SELECT MAF_tier, count(*) FROM variants GROUP BY MAF_tier ORDER BY MAF_tier"):
    print(f"  {tier}: {n:,}")

# 3) annotated file for inspection (a few key columns)
with gzip.open(OUT, "wt", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Variant_ID", "RSID", "Non_Finnish_European_AF", "Joint_AF", "MAF", "MAF_tier"])
    for vid, rs, nfe, jaf, tier in db.execute(
            "SELECT Variant_ID, RSID, Non_Finnish_European_AF, Joint_AF, MAF_tier FROM variants"):
        maf_v = None if nfe is None else round(min(nfe, 1 - nfe), 6)
        w.writerow([vid, rs, nfe, jaf, maf_v, tier])
import os
print(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.0f} MB); done in {time.time()-t0:.0f}s")
db.close()
