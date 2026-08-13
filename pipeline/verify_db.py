#!/usr/bin/env python3
"""Battery of sanity checks on evo2.db (June-18 gnomAD v4.1 joint schema)."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3

db = sqlite3.connect(str(config.DB_PATH))
db.row_factory = sqlite3.Row
q = lambda s: db.execute(s).fetchall()
def hr(t): print("\n" + "=" * 70 + f"\n{t}\n" + "=" * 70)

hr("ROW COUNT / UNIQUENESS")
print("total rows           :", q("SELECT count(*) c FROM variants")[0]["c"])
print("distinct Variant_ID  :", q("SELECT count(DISTINCT Variant_ID) c FROM variants")[0]["c"])
print("distinct RSID        :", q("SELECT count(DISTINCT RSID) c FROM variants")[0]["c"])

hr("CHROMOSOME COVERAGE (count per Chromosome)")
for r in q("SELECT Chromosome, count(*) n, min(Position) lo, max(Position) hi FROM variants GROUP BY Chromosome ORDER BY n DESC"):
    print(f"  chr {str(r['Chromosome']):<4} {r['n']:>9,}   pos {r['lo']:>12,} .. {r['hi']:>12,}")

hr("NULL HANDLING (key columns)")
for col in ["Evo2_7B_NoRC_Delta", "Evo2_40B_NoRC_Delta", "PhastCons_100way",
            "PhyloP_100way", "Joint_AF", "Non_Finnish_European_AF", "Gene", "Amino_Acid_Change"]:
    n = q(f'SELECT count(*) c FROM variants WHERE "{col}" IS NULL')[0]["c"]
    print(f"  NULL in {col:<26}: {n:>10,}")

hr("SCORE DISTRIBUTIONS (delta scores)")
for col in ["Evo2_7B_NoRC_Delta", "Evo2_40B_NoRC_Delta"]:
    r = q(f'SELECT min("{col}") lo, max("{col}") hi, avg("{col}") mu FROM variants')[0]
    print(f"  {col:<26}: min={r['lo']:.3f}  max={r['hi']:.3f}  mean={r['mu']:.4f}")

hr("ALLELE-FREQUENCY RANGES (gnomAD v4.1 joint)")
for col in ["Joint_AF", "Non_Finnish_European_AF", "East_Asian_AF", "African_AF", "South_Asian_AF"]:
    r = q(f'SELECT min("{col}") lo, max("{col}") hi, count("{col}") n FROM variants')[0]
    print(f"  {col:<26}: min={r['lo']}  max={r['hi']}  non-null={r['n']:,}")

for label, col in [("Functional_Annotation", "Functional_Annotation"),
                   ("Exonic_Function", "Exonic_Function"),
                   ("ENCODE_cCRE_Annotation", "ENCODE_cCRE_Annotation"),
                   ("Repeat_Class", "Repeat_Class")]:
    hr(f"VALUE COUNTS: {label}")
    for r in q(f'SELECT "{col}" v, count(*) n FROM variants GROUP BY "{col}" ORDER BY n DESC LIMIT 15'):
        print(f"  {str(r['v'])[:45]:<46} {r['n']:>10,}")

hr("TOP GENES BY VARIANT COUNT")
for r in q('SELECT Gene v, count(*) n FROM variants WHERE Gene IS NOT NULL GROUP BY Gene ORDER BY n DESC LIMIT 12'):
    print(f"  {str(r['v'])[:40]:<41} {r['n']:>9,}")

hr("FTS5 SEARCH TEST  (MATCH 'APOE')")
rows = q("SELECT v.RSID, v.Chromosome, v.Position, v.Gene, v.Evo2_7B_NoRC_Delta d "
         "FROM variants v JOIN variants_fts f ON v.rowid=f.rowid "
         "WHERE variants_fts MATCH 'APOE' LIMIT 5")
print("  matched (showing 5 of) ...")
for r in rows:
    print(f"    {r['RSID']}  chr{r['Chromosome']}:{r['Position']}  {r['Gene']}  7b_delta={r['d']}")

hr("POINT LOOKUP  (Variant_ID = '19:44908684:T:C'  = rs429358, APOE e4)")
for r in q("SELECT * FROM variants WHERE Variant_ID='19:44908684:T:C'"):
    d = dict(r)
    for k in ["Variant_ID", "RSID", "Joint_AF", "Non_Finnish_European_AF", "East_Asian_AF",
              "Evo2_7B_NoRC_Delta", "Evo2_40B_NoRC_Delta",
              "Functional_Annotation", "Gene", "ENCODE_cCRE_Annotation", "PhyloP_100way"]:
        print(f"    {k:<26}: {d[k]}")

hr("RANGE QUERY PLAN (uses idx_chr_pos?)")
for r in q("EXPLAIN QUERY PLAN SELECT * FROM variants WHERE Chromosome='1' AND Position BETWEEN 60000 AND 70000"):
    print("   ", r["detail"] if "detail" in r.keys() else tuple(r))

db.close()
