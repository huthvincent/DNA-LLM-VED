#!/usr/bin/env python3
"""Re-ingest GWAS associations from the pre-filtered LDSC_P_1E-4 tree.

The user supplied per-trait files already filtered to P<=1e-4 (format:
'rsID:Ref:Alt<TAB>P', no header) for essentially ALL LDSC traits — far more than
the streamed-zip download recovered. This replaces the gwas_studies / gwas_assoc /
gwas_import_log tables with the complete set, bridging /gwas (LDSC) and /coloc.
Variants are linked to our DB by rsID.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, os, csv, glob, re, time, math

DB = str(config.DB_PATH)
ROOT = str(config.GWAS_P_DIR)
ANNO = str(config.GWAS_ANNO_DIR)
LDSC_JSON = str(config.GWAS_LDSC_JSON)


def load_tsv(p):
    with open(p) as fh:
        r = csv.reader(fh, delimiter="\t"); next(r); return [row for row in r]


def num(x):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return None


import json
ldsc = json.load(open(LDSC_JSON))["sources"]
ukb_cat = {r["id"]: r["category"] for r in ldsc.get("UKB", [])}
mvp = {r[0]: {"name": r[5], "cat": r[4], "n": num(r[6]), "cases": num(r[7]), "controls": num(r[8])} for r in load_tsv(f"{ANNO}/MVP_annotation.txt") if r}
ukb = {r[1]: {"name": r[2], "tag": r[0]} for r in load_tsv(f"{ANNO}/UKB_annotation.txt") if len(r) > 3}
finn = {r[1]: {"name": r[0], "cat": r[3], "cases": num(r[4]), "controls": num(r[6])} for r in load_tsv(f"{ANNO}/FinnGen_annotation.txt") if len(r) > 6}
mvp_keys = set(mvp)
print(f"annotations: MVP {len(mvp)}, UKB {len(ukb)}, FinnGen {len(finn)}")


def ancestry_of(base):
    for a in ("EUR", "AFR", "AMR", "EAS", "SAS", "ASJ", "FIN", "MID", "META", "AGR"):
        if f".{a}." in base or f"_{a}_" in base or f"_{a}." in base or base.endswith(f"_{a}"):
            return a
    return None


def meta_for(cohort, base):
    anc = ancestry_of(base)
    stem = base
    for suf in ("_1E-4.txt", ".tsv.gz", ".txt.gz", ".gz", ".txt"):
        if stem.endswith(suf):
            stem = stem[:-len(suf)]; break
    if cohort == "MVP":
        toks = stem.split(".")
        code = next((t for t in toks if t in mvp_keys), toks[2] if len(toks) > 2 else stem)
        a = mvp.get(code)
        if a:
            return a["name"], a["cat"], a["n"], a["cases"], a["controls"], anc, code
        return code, None, None, None, None, anc, None
    if cohort == "FinnGen":
        code = re.sub(r"^finngen_R\d+_", "", stem)
        a = finn.get(code)
        if a:
            n = (a["cases"] or 0) + (a["controls"] or 0) or None
            return a["name"], a["cat"], n, a["cases"], a["controls"], anc or "FIN", code
        return code, None, None, None, None, anc or "FIN", None
    if cohort == "UKB":
        acc = stem
        a = ukb.get(acc)
        if a:
            return a["name"], ukb_cat.get(acc), None, None, None, ancestry_of(a["tag"]) or anc, acc
        return acc, ukb_cat.get(acc), None, None, None, anc, acc
    return stem, None, None, None, None, anc, None


t0 = time.time()
db = sqlite3.connect(DB)
db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA cache_size=-400000;")
print("loading rsID set ...")
rsids = set(r[0] for r in db.execute("SELECT RSID FROM variants WHERE RSID IS NOT NULL"))
print(f"  {len(rsids):,} rsIDs ({time.time()-t0:.0f}s)")

db.executescript("""
DROP TABLE IF EXISTS gwas_studies; DROP TABLE IF EXISTS gwas_assoc; DROP TABLE IF EXISTS gwas_import_log;
CREATE TABLE gwas_studies (study_id INTEGER PRIMARY KEY, cohort TEXT, trait_name TEXT, trait_file_name TEXT,
  phenotype_category TEXT, description TEXT, sample_size INTEGER, case_number INTEGER, control_number INTEGER,
  ancestry TEXT, pmid TEXT, publication TEXT, publication_year INTEGER, source_zip_file TEXT, notes TEXT);
CREATE TABLE gwas_assoc (study_id INTEGER, rsid TEXT, p_value REAL);
CREATE TABLE gwas_import_log (study_id INTEGER, source_file TEXT, total_variants INTEGER, matched_variants INTEGER,
  unmatched_variants INTEGER, matched_rate REAL, min_p REAL, n_p_lt_5e_8 INTEGER, n_p_lt_1e_5 INTEGER,
  import_time REAL, status TEXT, message TEXT);
""")

study_id = 0
for cohort in ("UKB", "FinnGen", "MVP"):
    files = sorted(glob.glob(f"{ROOT}/{cohort}/*.txt"))
    t_co = time.time()
    for path in files:
        base = os.path.basename(path)
        total = 0; matched = []; minp = 1.0; n8 = 0
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                total += 1
                try:
                    pv = float(parts[1])
                except ValueError:
                    continue
                if pv < minp:
                    minp = pv
                if pv < 5e-8:
                    n8 += 1
                rs = parts[0].split(":", 1)[0]
                if rs in rsids:
                    matched.append((study_id, rs, pv))
        db.executemany("INSERT INTO gwas_assoc VALUES (?,?,?)", matched)
        name, cat, n, cases, controls, anc, note = meta_for(cohort, base)
        db.execute("INSERT INTO gwas_studies (study_id,cohort,trait_name,trait_file_name,phenotype_category,"
                   "description,sample_size,case_number,control_number,ancestry,source_zip_file,notes) "
                   "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                   (study_id, cohort, name, base, cat, name, n, cases, controls, anc, f"LDSC_P_1E-4/{cohort}", note))
        rate = len(matched) / total if total else None
        db.execute("INSERT INTO gwas_import_log VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                   (study_id, base, total, len(matched), total - len(matched), rate, minp, n8, total,
                    time.time(), "ok", None))
        study_id += 1
    db.commit()
    print(f"  {cohort}: {len(files)} traits  ({time.time()-t_co:.0f}s)", flush=True)

print("[index] ...", flush=True)
db.execute("CREATE INDEX idx_assoc_study ON gwas_assoc(study_id)")
db.execute("CREATE INDEX idx_assoc_rsid ON gwas_assoc(rsid)")
db.execute("CREATE INDEX idx_studies_cohort ON gwas_studies(cohort, phenotype_category)")
db.commit()
ns = db.execute("SELECT count(*) FROM gwas_studies").fetchone()[0]
na = db.execute("SELECT count(*) FROM gwas_assoc").fetchone()[0]
print(f"[done] {ns:,} studies, {na:,} associations in {time.time()-t0:.0f}s; db={os.path.getsize(DB)/1e9:.2f}GB")
db.close()
