#!/usr/bin/env python3
"""Re-annotate GWAS trait names from the authoritative annotation files (July-16 fix).

Bug: the S-LDSC table (gwas_ldsc.json, built by build_gwas.py from raw LDSC files by
COLUMN POSITION) picked the wrong column for FinnGen trait names -- 2,042 of 2,297
FinnGen traits carried a name belonging to a different endpoint (e.g. every Alzheimer
endpoint collapsed to "Alzheimer's disease"; ANTIDEPRESSANTS -> "Injury, NOS").

The z / coef / se / p / n statistics are keyed by trait id and are unaffected (MVP
validates 100%), so we re-map only the display name + category, by id, from:

  FinnGen : FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx   NAME -> LONGNAME
            (the .txt supplied under that name is byte-identical to the MVP file,
             so the xlsx is the real FinnGen source; the R12 annotation we already
             ship is used as a fallback for the ~65 endpoints absent from the xlsx,
             and supplies the phenotype category)
  MVP     : MVP_GWAS_annotation.txt                                Trait -> Description
  UKB WGS : UKB_WGS_GWAS_annotation.tsv                            accessionId -> reportedTrait

The same mapping is applied to gwas_studies.trait_name (the /coloc + trait tables) so
every module displays identical names for the same trait id.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import json, csv, sqlite3, shutil, sys, os
import openpyxl

FIX = str(config.GWAS_ANNO_DIR)
DB = str(config.DB_PATH)
LDSC_JSON = str(config.GWAS_LDSC_JSON)
R12 = str(config.GWAS_ANNO_DIR / "FinnGen_annotation.txt")
APPLY = "--apply" in sys.argv


# ---------------------------------------------------------------- authoritative maps
def finngen_map():
    """NAME -> (LONGNAME, category). xlsx is authoritative for the name; the R12
    annotation fills endpoints missing from the xlsx and provides the category."""
    ws = openpyxl.load_workbook(f"{FIX}/FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx",
                                read_only=True)["Sheet 1"]
    rows = ws.iter_rows(values_only=True)
    H = {h: i for i, h in enumerate(next(rows)) if h}
    long_ = {}
    for r in rows:
        nm, ln = r[H["NAME"]], r[H["LONGNAME"]]
        if nm and not str(nm).startswith("#") and ln:
            long_[str(nm).strip()] = str(ln).strip()
    r12 = {}
    with open(R12) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r12[r["phenocode"].strip()] = (r["phenotype"].strip(), r["category"].strip())
    out = {}
    for code in set(long_) | set(r12):
        name = long_.get(code) or (r12[code][0] if code in r12 else None)
        cat = r12[code][1] if code in r12 else None
        if name:
            out[code] = (name, cat)
    return out, len(long_), len(r12)


def mvp_map():
    out = {}
    with open(f"{FIX}/MVP_GWAS_annotation.txt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = (r["Trait"] or "").strip()
            if t:
                out[t] = ((r["Description"] or "").strip() or t, (r["SubCategory"] or "").strip() or None)
    return out


def ukb_map():
    out = {}
    with open(f"{FIX}/UKB_WGS_GWAS_annotation.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a = (r.get("accessionId") or "").strip()
            nm = (r.get("reportedTrait") or "").strip()
            if a and nm:
                out[a] = (nm, None)      # keep the existing phecode category for UKB
    return out


fin, n_xlsx, n_r12 = finngen_map()
mvp, ukb = mvp_map(), ukb_map()
REF = {"FinnGen": fin, "MVP": mvp, "UKB": ukb}
print(f"authoritative: FinnGen {len(fin)} (xlsx {n_xlsx} + R12 {n_r12}), MVP {len(mvp)}, UKB {len(ukb)}\n")

report = {}

# ---------------------------------------------------------------- 1) gwas_ldsc.json
d = json.load(open(LDSC_JSON))
for co, recs in d["sources"].items():
    ref = REF[co]
    fixed = cat_fixed = unmapped = 0
    unmapped_ids = []
    for r in recs:
        key = str(r.get("id", "")).strip()
        if key not in ref:
            unmapped += 1
            if len(unmapped_ids) < 20:
                unmapped_ids.append(key)
            continue
        name, cat = ref[key]
        if name and name != r.get("name"):
            r["name"] = name; fixed += 1
        if cat and cat != r.get("category"):
            r["category"] = cat; cat_fixed += 1
    report[co] = {"total": len(recs), "name_corrected": fixed, "category_corrected": cat_fixed,
                  "unmapped": unmapped, "unmapped_ids": unmapped_ids}
    print(f"[ldsc {co}] total={len(recs)} name_corrected={fixed} category_corrected={cat_fixed} unmapped={unmapped}")
d["version"] = "2026-07-24-annotation-fix"
d["annotation_sources"] = {
    "FinnGen": "FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx (NAME->LONGNAME); R12 annotation fallback + category",
    "MVP": "MVP_GWAS_annotation.txt (Trait->Description)",
    "UKB": "UKB_WGS_GWAS_annotation.tsv (accessionId->reportedTrait)"}

# ---------------------------------------------------------------- 2) gwas_studies
db = sqlite3.connect(DB)
studies = {}
for co in REF:
    rows = db.execute("SELECT study_id, notes, trait_name FROM gwas_studies WHERE cohort=?", (co,)).fetchall()
    ref = REF[co]
    upd = []
    unmapped = 0
    for sid, code, cur in rows:
        key = (code or "").strip()
        a = ref.get(key)
        if not a:
            unmapped += 1
            continue
        if a[0] and a[0] != (cur or ""):
            upd.append((a[0], sid))
    studies[co] = {"total": len(rows), "name_corrected": len(upd), "unmapped": unmapped}
    print(f"[studies {co}] total={len(rows)} name_corrected={len(upd)} unmapped={unmapped}")
    if APPLY and upd:
        db.executemany("UPDATE gwas_studies SET trait_name=? WHERE study_id=?", upd)
report["gwas_studies"] = studies

if APPLY:
    shutil.copy(LDSC_JSON, LDSC_JSON + ".bak-preannofix")
    json.dump(d, open(LDSC_JSON, "w"))
    db.commit()
    json.dump(report, open(f"{BASE}/evo2_database/gwas_annotation_fix_report.json", "w"), indent=1)
    print("\n[applied] gwas_ldsc.json rewritten (.bak-preannofix kept), gwas_studies updated, report written")
else:
    print("\n[dry-run] nothing written -- rerun with --apply")
db.close()
