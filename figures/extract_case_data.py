#!/usr/bin/env python3
"""Pull region-scoped aggregates/samples from evo2.db for the NAR Fig 3 case study.

All four panels (a-d) are drawn from ONE region: the APOE locus +/- 1 Mb on chr19.
Because the database holds only common (NFE MAF >= 5%) variants, protein-truncating
consequences (stopgain/stoploss/startloss) never occur in any single locus (LoF
cannot reach 5% frequency), so coding consequence is limited to synonymous /
nonsynonymous. We keep whatever annotation classes are actually present in the
region (with a minimum n for the boxplot panels) and drop the empty/trivial ones.

Schema = the June-18 gnomAD v4.1 joint re-annotation (new human-readable names).
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, json

DB = str(config.DB_PATH)
REGION = {"chr": "19", "start": 43900000, "end": 45900000}  # APOE locus +/- 1 Mb (~4.8k common variants)
CORE   = {"start": 44880000, "end": 44950000}               # APOE 167-site core (highlighted in panel d)
MIN_N  = 10                                                  # drop annotation classes thinner than this (panels a/b)

db = sqlite3.connect(DB); db.row_factory = sqlite3.Row
W = (REGION["chr"], REGION["start"], REGION["end"])
q = lambda s, p=(): [dict(r) for r in db.execute(s, p).fetchall()]
out = {"region": REGION, "core": CORE}

REG = "Chromosome=? AND Position BETWEEN ? AND ? AND Evo2_40B_NoRC_Delta IS NOT NULL"

# A. delta by functional region (Functional_Annotation) -- keep present classes with n >= MIN_N
forder = ["exonic", "UTR5", "UTR3", "ncRNA_exonic", "upstream", "downstream",
          "intronic", "ncRNA_intronic", "intergenic"]
out["delta_by_func"] = {}
for fc in forder:
    rows = q(f"SELECT Evo2_40B_NoRC_Delta d FROM variants WHERE {REG} AND Functional_Annotation=?", (*W, fc))
    if len(rows) >= MIN_N:
        out["delta_by_func"][fc] = [r["d"] for r in rows]

# B. delta by coding consequence (Exonic_Function) -- keep only classes actually present in the region
eorder = ["synonymous SNV", "nonsynonymous SNV", "startloss", "stoploss", "stopgain",
          "nonframeshift substitution", "frameshift substitution"]
out["delta_by_exonic"] = {}
for ef in eorder:
    rows = q(f"SELECT Evo2_40B_NoRC_Delta d FROM variants WHERE {REG} AND Exonic_Function=?", (*W, ef))
    if rows:
        out["delta_by_exonic"][ef] = [r["d"] for r in rows]

# C. all 6 model x strand-strategy delta scores for region variants (Spearman computed in the fig script)
dcols = ["Evo2_7B_NoRC_Delta", "Evo2_7B_AvgRC_Delta", "Evo2_7B_WeightedRC_Delta",
         "Evo2_40B_NoRC_Delta", "Evo2_40B_AvgRC_Delta", "Evo2_40B_WeightedRC_Delta"]
out["region_scores"] = q(
    f"SELECT {','.join(dcols)} FROM variants WHERE Chromosome=? AND Position BETWEEN ? AND ? AND "
    + " AND ".join(f"{c} IS NOT NULL" for c in dcols), W)

# D. region variant landscape (all scored common variants in the +/-1 Mb region).
#    Aliased back to CHR/BP/d7/d40/ef/fc so the plotting code is schema-agnostic.
out["gene_case"] = {"gene": "APOE locus", "region": REGION, "core": CORE, "rows": q(
    "SELECT Chromosome AS CHR, Position AS BP, Evo2_7B_NoRC_Delta AS d7, Evo2_40B_NoRC_Delta AS d40, "
    "Exonic_Function AS ef, Functional_Annotation AS fc FROM variants "
    f"WHERE {REG} ORDER BY Position", W)}

# context: functional-class counts within the region
out["func_counts"] = q(
    "SELECT Functional_Annotation f, count(*) n FROM variants WHERE Chromosome=? AND Position BETWEEN ? AND ? "
    "GROUP BY f ORDER BY n DESC", W)

config.CASE_DATA_JSON.parent.mkdir(parents=True, exist_ok=True)
json.dump(out, open(str(config.CASE_DATA_JSON), "w"))
print("region:", REGION, "core:", CORE)
print("A delta_by_func:", {k: len(v) for k, v in out["delta_by_func"].items()})
print("B delta_by_exonic:", {k: len(v) for k, v in out["delta_by_exonic"].items()})
print("C region_scores n:", len(out["region_scores"]))
print("D landscape n:", len(out["gene_case"]["rows"]))
db.close()
