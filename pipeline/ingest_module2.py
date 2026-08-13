#!/usr/bin/env python3
"""Ingest Module-2 auxiliary data into evo2.db:

  * allele_age  -- GEVA (AgeMean_Jnt) + SDS per variant, keyed by RSID.
                   NOTE: the GEVA/SDS source files are GRCh37/hg19, so their
                   coordinates do NOT match our GRCh38 db; the rsID matches 100%,
                   so we join by rsID. We keep only rsIDs present in our db.
  * region_sets -- built-in genomic region sets (currently HAR, hg38 BED) used by
                   the Region-Set analysis (intersected with variants by position).

Adds two small tables; the 6.48M-row `variants` table is left untouched.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import gzip, sqlite3, glob, time, os

DB  = str(config.DB_PATH)
AA  = str(config.SOURCE_ROOT / "allele_age")
HAR = str(config.HAR_BED)

t0 = time.time()
db = sqlite3.connect(DB)
db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA cache_size=-400000;")

print("[main] loading db rsIDs into a set ...")
main_rsids = set(r[0] for r in db.execute("SELECT RSID FROM variants WHERE RSID IS NOT NULL"))
print(f"[main] {len(main_rsids):,} rsIDs")

db.execute("DROP TABLE IF EXISTS allele_age")
db.execute("CREATE TABLE allele_age (rsid TEXT PRIMARY KEY, geva REAL, sds REAL)")


def stream(path, valcol, set_col):
    f = gzip.open(path, "rt")
    h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
    ri, vi = ix["rsid"], ix[valcol]
    buf, n, kept = [], 0, 0
    sql = (f"INSERT INTO allele_age(rsid,{set_col}) VALUES(?,?) "
           f"ON CONFLICT(rsid) DO UPDATE SET {set_col}=excluded.{set_col}")
    for line in f:
        r = line.rstrip("\n").split("\t")
        rs = r[ri]
        if not rs or rs == "NA" or rs not in main_rsids:
            continue
        try:
            val = float(r[vi])
        except (ValueError, IndexError):
            continue
        buf.append((rs, val)); kept += 1
        if len(buf) >= 50000:
            db.executemany(sql, buf); buf.clear()
        n += 1
    if buf:
        db.executemany(sql, buf)
    db.commit()
    return kept


geva_file = glob.glob(AA + "/*GEVA*.txt.gz")[0]
sds_file  = glob.glob(AA + "/*SDS*.txt.gz")[0]
print("[geva] ingesting AgeMean_Jnt ...")
g = stream(geva_file, "AgeMean_Jnt", "geva"); print(f"[geva] kept {g:,}  ({time.time()-t0:.0f}s)")
print("[sds]  ingesting SDS ...")
s = stream(sds_file, "SDS", "sds"); print(f"[sds]  kept {s:,}  ({time.time()-t0:.0f}s)")

# region sets (HAR, hg38)
db.execute("DROP TABLE IF EXISTS region_sets")
db.execute("CREATE TABLE region_sets (set_name TEXT, chrom TEXT, start INTEGER, end INTEGER)")
har = []
for line in open(HAR):
    p = line.split()
    if len(p) < 3:
        continue
    c = p[0].replace("chr", "")
    try:
        har.append(("HAR", c, int(p[1]), int(p[2])))
    except ValueError:
        continue
db.executemany("INSERT INTO region_sets VALUES (?,?,?,?)", har)
db.execute("CREATE INDEX idx_regionset ON region_sets(set_name, chrom, start, end)")
db.commit()
print(f"[region_sets] HAR: {len(har):,} regions")

# report join coverage
n_aa = db.execute("SELECT count(*) FROM allele_age").fetchone()[0]
n_geva = db.execute("SELECT count(*) FROM allele_age WHERE geva IS NOT NULL").fetchone()[0]
n_sds = db.execute("SELECT count(*) FROM allele_age WHERE sds IS NOT NULL").fetchone()[0]
cov = db.execute("SELECT count(*) FROM variants v JOIN allele_age a ON a.rsid=v.RSID").fetchone()[0]
print(f"[allele_age] rows={n_aa:,}  geva={n_geva:,}  sds={n_sds:,}  | variants with any allele-age: {cov:,}")
print("[vacuum] ...")
db.execute("VACUUM"); db.close()
print(f"[done] {time.time()-t0:.0f}s; db size = {os.path.getsize(DB)/1e9:.2f} GB")
