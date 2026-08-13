#!/usr/bin/env python3
"""Rebuild the allele_age table from the July-12 corrected GEVA + SDS files, and add
the AF_bin column used by the updated (subsampling) allele-age analyses.

Why this replaces ingest_module2.py's allele_age:
  * both new files are derived from the SAME source as the database itself, so the
    allele-age annotation is now consistent with the variant table (previously the
    GEVA/SDS files were hg19-derived and joined by rsID from a different release);
  * the SDS file now also contains variants whose REF/ALT were reversed during allele
    harmonisation (record count roughly doubled), and carries `SDS_Final` -- the SDS
    value re-signed to match the allele orientation used in this database. We store
    `SDS_Final` as `sds`; the raw unharmonised SDS is deliberately NOT exposed.
  * GEVA supplies both AgeMean_Jnt and AgeMedian_Jnt. Per the maintainer's decision the
    database standardises on **AgeMean_Jnt** (stored as `geva`) for every analysis;
    AgeMedian_Jnt is kept alongside as `age_median` for reference/validation.

AF_bin = equal-width bins of Non_Finnish_European_AF (verified identical to the `AF_bin`
column shipped in both files): AF_0-0.2 / AF_0.2-0.4 / AF_0.4-0.6 / AF_0.6-0.8 / AF_0.8-1.0.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, gzip, time, os, sys

DB = str(config.DB_PATH)
GEVA = str(config.GEVA_TSV)
SDS = str(config.SDS_TSV)
AF_LEVELS = ["AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0"]


def scan(path, want):
    """Stream a .gz TSV, yielding the requested columns by header name."""
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = [hdr.index(c) for c in want]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= ix[-1]:
                continue
            yield [f[i] for i in ix]


def fnum(x):
    try:
        v = float(x)
        return v if v == v else None
    except (ValueError, TypeError):
        return None


t0 = time.time()

# ---------------------------------------------------------------- GEVA
age = {}
for rs, mean, med in scan(GEVA, ["RSID", "AgeMean_Jnt", "AgeMedian_Jnt"]):
    if rs and rs != "." and rs not in age:
        age[rs] = (fnum(mean), fnum(med))
print(f"GEVA rows keyed by rsID: {len(age):,}  ({time.time()-t0:.0f}s)", flush=True)

# ---------------------------------------------------------------- SDS (SDS_Final)
sds = {}
flips = {"Match": 0, "other": 0}
for rs, final, match in scan(SDS, ["RSID", "SDS_Final", "SDS_allele_match"]):
    flips["Match" if match == "Match" else "other"] = flips.get("Match" if match == "Match" else "other", 0) + 1
    if rs and rs != "." and rs not in sds:
        v = fnum(final)
        if v is not None:
            sds[rs] = v
print(f"SDS rows keyed by rsID: {len(sds):,}  (allele_match: {flips})  ({time.time()-t0:.0f}s)", flush=True)

# ---------------------------------------------------------------- write allele_age
db = sqlite3.connect(DB)
db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA cache_size=-400000;")
db.executescript("""
DROP TABLE IF EXISTS allele_age;
CREATE TABLE allele_age (rsid TEXT PRIMARY KEY, geva REAL, age_median REAL, sds REAL);
""")
rows = []
for rs in set(age) | set(sds):
    m, med = age.get(rs, (None, None))
    rows.append((rs, m, med, sds.get(rs)))
db.executemany("INSERT INTO allele_age VALUES (?,?,?,?)", rows)
db.commit()
n_g = db.execute("SELECT count(*) FROM allele_age WHERE geva IS NOT NULL").fetchone()[0]
n_s = db.execute("SELECT count(*) FROM allele_age WHERE sds IS NOT NULL").fetchone()[0]
print(f"allele_age: {len(rows):,} rsIDs  (geva {n_g:,}, sds {n_s:,})  ({time.time()-t0:.0f}s)", flush=True)

# how many actually join to our variant table
j_g, j_s = db.execute("""SELECT sum(a.geva IS NOT NULL), sum(a.sds IS NOT NULL)
                         FROM variants v JOIN allele_age a ON a.rsid=v.RSID""").fetchone()
print(f"  joined to variants: geva {j_g:,}, sds {j_s:,}")

# ---------------------------------------------------------------- AF_bin on variants
cols = [r[1] for r in db.execute("PRAGMA table_info(variants)")]
if "AF_bin" not in cols:
    db.execute("ALTER TABLE variants ADD COLUMN AF_bin TEXT")
af = "Non_Finnish_European_AF"
db.execute(f"""
UPDATE variants SET AF_bin = CASE
  WHEN {af} IS NULL   THEN NULL
  WHEN {af} <  0.2    THEN '{AF_LEVELS[0]}'
  WHEN {af} <  0.4    THEN '{AF_LEVELS[1]}'
  WHEN {af} <  0.6    THEN '{AF_LEVELS[2]}'
  WHEN {af} <  0.8    THEN '{AF_LEVELS[3]}'
  ELSE '{AF_LEVELS[4]}' END
""")
db.execute("CREATE INDEX IF NOT EXISTS idx_af_bin ON variants(AF_bin)")
db.commit()
print("AF_bin counts:")
for b, n in db.execute("SELECT AF_bin, count(*) FROM variants GROUP BY AF_bin ORDER BY AF_bin"):
    print(f"  {b}: {n:,}")
print(f"[done] {time.time()-t0:.0f}s; db={os.path.getsize(DB)/1e9:.2f}GB")
db.close()
