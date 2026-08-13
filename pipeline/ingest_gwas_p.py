#!/usr/bin/env python3
"""Ingest raw GWAS P-values (P <= 1e-4) into evo2.db for the co-localization module.

For each cohort zip (streamed member-by-member with libarchive; members are
gzipped 'rsID:Ref:Alt<TAB>P' files), awk does the heavy per-line work (count,
min_p, n<5e-8, n<1e-5, emit P<=1e-4 rows + a 20k-row rsID sample) so Python only
handles the small kept/sample sets. Variants are linked by rsID (the GWAS files
have no coordinates; our DB has unique rsID). Builds three tables:
  gwas_studies (metadata), gwas_assoc (study_id, rsid, p_value), gwas_import_log.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import sqlite3, subprocess, gzip, time, os, csv, glob, re, zlib, struct


def iter_members(path, bufsize=1 << 24):
    """Yield (filename, decompressed_bytes) for each STORED gzip member of a
    streamed/no-central-directory zip. Error-placeholder (.gz_Error.txt) and
    other non-.gz members yield (filename, None). Bypasses libarchive, which
    chokes on the mixed STORED-gzip + placeholder layout these zips use."""
    SIG = b"PK\x03\x04"
    f = open(path, "rb"); buf = f.read(bufsize)
    while True:
        i = buf.find(SIG)
        while i < 0:
            m = f.read(bufsize)
            if not m:
                return
            buf = buf[-3:] + m; i = buf.find(SIG)
        buf = buf[i:]
        while len(buf) < 30:
            m = f.read(bufsize)
            if not m:
                return
            buf += m
        fnl, enl = struct.unpack("<HH", buf[26:30])
        while len(buf) < 30 + fnl + enl:
            m = f.read(bufsize)
            if not m:
                return
            buf += m
        fname = buf[30:30 + fnl].decode("utf-8", "replace"); body = buf[30 + fnl + enl:]
        if fname.endswith(".gz"):
            d = zlib.decompressobj(32 + 15); out = bytearray(); cur = body; ok = True
            while True:
                if not cur:
                    cur = f.read(bufsize)
                    if not cur:
                        ok = False; break
                try:
                    out += d.decompress(cur)
                except zlib.error:
                    ok = False; break
                if d.eof:
                    buf = d.unused_data; break
                cur = b""
            if ok and d.eof:
                yield fname, bytes(out); continue
            buf = cur if cur else b""
        else:
            yield fname, None; buf = body

DB = str(config.DB_PATH)
GZDIR = str(config.SOURCE_ROOT / "gwas" / "GWAS_P")
ANNO = str(config.GWAS_ANNO_DIR)
LDSC_JSON = str(config.GWAS_LDSC_JSON)
P_KEEP = 1e-4
SAMPLE = 20000

AWK = r'''BEGIN{tot=0; n8=0; n4=0; minp=1}
{ tot++; p=$2+0; if(tot==1||p<minp)minp=p; if(p<5e-8)n8++;
  if(p<1e-4){n4++; print "K\t" $1 "\t" $2}
  if(tot<=20000){r=$1; sub(/:.*/,"",r); print "S\t" r} }
END{ print "T\t" tot "\t" n8 "\t" n4 "\t" minp }'''


def load_tsv(path):
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t"); hdr = next(r)
        return hdr, [row for row in r]


def num(x):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return None


# ---- annotation lookups -----------------------------------------------------
import json
ldsc = json.load(open(LDSC_JSON))["sources"]
ukb_cat = {r["id"]: r["category"] for r in ldsc.get("UKB", [])}   # by Study.Accession

# MVP: keyed by Trait(0); name=Description(5), cat=SubCategory(4), n=num_samples.META(6), cases(7), controls(8)
_h, mvp_rows = load_tsv(f"{ANNO}/MVP_annotation.txt")
mvp = {r[0]: {"name": r[5], "cat": r[4], "n": num(r[6]), "cases": num(r[7]), "controls": num(r[8])} for r in mvp_rows if r}
# UKB: files are named <Study Accession>.tsv.gz; key by accession(1). name=Reported trait(2), tag(0)
_h, ukb_rows = load_tsv(f"{ANNO}/UKB_annotation.txt")
ukb = {r[1]: {"name": r[2], "tag": r[0], "acc": r[1]} for r in ukb_rows if len(r) > 3}
# FinnGen: keyed by phenocode(1); name=phenotype(0), cat=category(3), cases(4), controls(6)
_h, finn_rows = load_tsv(f"{ANNO}/FinnGen_annotation.txt")
finn = {r[1]: {"name": r[0], "cat": r[3], "cases": num(r[4]), "controls": num(r[6])} for r in finn_rows if len(r) > 6}
mvp_keys = set(mvp); finn_keys = set(finn)
print(f"annotations: MVP {len(mvp)}, UKB {len(ukb)}, FinnGen {len(finn)}")


def ancestry_of(base):
    for a in ("EUR", "AFR", "AMR", "EAS", "SAS", "ASJ", "FIN", "MID", "META", "AGR"):
        if f".{a}." in base or f"_{a}_" in base or f"_{a}." in base or base.endswith(f"_{a}"):
            return a
    return None


def meta_for(cohort, base):
    """Return (trait_name, category, sample_size, cases, controls, ancestry, notes)."""
    anc = ancestry_of(base)
    stem = base
    for suf in (".tsv.gz", ".txt.gz", ".gz"):
        if stem.endswith(suf):
            stem = stem[:-len(suf)]
    if cohort == "MVP":
        toks = stem.split(".")
        code = next((t for t in toks if t in mvp_keys), toks[2] if len(toks) > 2 else stem)
        a = mvp.get(code)
        if a:
            return a["name"], a["cat"], a["n"], a["cases"], a["controls"], anc, code
        return code, None, None, None, None, anc, None
    if cohort == "FinnGen":
        code = re.sub(r"^finngen_R\d+_", "", stem)   # filename = finngen_R12_<PHENOCODE>
        a = finn.get(code)
        if a:
            n = (a["cases"] or 0) + (a["controls"] or 0) or None
            return a["name"], a["cat"], n, a["cases"], a["controls"], anc or "FIN", code
        return code, None, None, None, None, anc or "FIN", None
    if cohort == "UKB":
        acc = stem  # filename = <Study Accession>.tsv.gz
        a = ukb.get(acc)
        if a:
            return a["name"], ukb_cat.get(acc), None, None, None, ancestry_of(a["tag"]) or anc, acc
        return acc, ukb_cat.get(acc), None, None, None, anc, acc
    return stem, None, None, None, None, anc, None


# ---- DB setup ---------------------------------------------------------------
db = sqlite3.connect(DB)
db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA cache_size=-400000;")
print("loading rsID set ...")
t0 = time.time()
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

ZIPS = [("MVP", f"{GZDIR}/MVP.zip"), ("FinnGen", f"{GZDIR}/FinnGen.zip"), ("UKB", f"{GZDIR}/UKB.zip")]
study_id = 0
t_all = time.time()
for cohort, zippath in ZIPS:
    if not os.path.exists(zippath):
        print(f"[skip] {zippath} not present"); continue
    print(f"\n=== {cohort}  ({os.path.getsize(zippath)/1e9:.1f}GB) ===", flush=True)
    nfiles = 0; nerr = 0; t_co = time.time()
    if True:
        for fname, data in iter_members(zippath):
            base = os.path.basename(fname)
            if data is None:
                db.execute("INSERT INTO gwas_import_log(study_id,source_file,status,message) VALUES(?,?,?,?)",
                           (study_id, base, "no_data", "download-error placeholder"))
                study_id += 1; nerr += 1; continue
            txt = data
            out = subprocess.run(["awk", "-F", "\t", AWK], input=txt, capture_output=True).stdout
            kept = []; samp = []; stats = None
            for ln in out.split(b"\n"):
                if ln[:2] == b"K\t":
                    parts = ln.split(b"\t"); kept.append((parts[1], parts[2]))
                elif ln[:2] == b"S\t":
                    samp.append(ln[2:])
                elif ln[:2] == b"T\t":
                    stats = ln.split(b"\t")
            try:
                total = int(stats[1]); n8 = int(stats[2]); n4 = int(stats[3]); minp = float(stats[4])
            except (ValueError, TypeError, IndexError, AttributeError):
                db.execute("INSERT INTO gwas_import_log(study_id,source_file,status,message) VALUES(?,?,?,?)",
                           (study_id, base, "skip", f"unparseable stats {stats!r}"[:150]))
                study_id += 1; continue
            ins = []
            for idv, p in kept:
                rs = idv.split(b":", 1)[0].decode()
                if rs in rsids:
                    ins.append((study_id, rs, float(p)))
            sm = sum(1 for r in samp if r.decode() in rsids)
            rate = sm / len(samp) if samp else None
            matched_est = int(round(total * rate)) if rate is not None else None
            db.executemany("INSERT INTO gwas_assoc VALUES (?,?,?)", ins)
            name, cat, n, cases, controls, anc, note = meta_for(cohort, base)
            db.execute("INSERT INTO gwas_studies (study_id,cohort,trait_name,trait_file_name,phenotype_category,"
                       "description,sample_size,case_number,control_number,ancestry,source_zip_file,notes) "
                       "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                       (study_id, cohort, name, base, cat, name, n, cases, controls, anc, f"{cohort}.zip", note))
            db.execute("INSERT INTO gwas_import_log VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                       (study_id, base, total, matched_est, (total - matched_est) if matched_est is not None else None,
                        rate, minp, n8, n4, time.time(), "ok", f"kept_imported={len(ins)}"))
            study_id += 1; nfiles += 1
            if nfiles % 200 == 0:
                db.commit()
                print(f"  {cohort}: {nfiles} traits  ({time.time()-t_co:.0f}s)  last='{name[:30]}' kept={len(ins)}", flush=True)
    db.commit()
    print(f"  {cohort} DONE: {nfiles} traits ({nerr} download-error placeholders) in {time.time()-t_co:.0f}s", flush=True)

print("\n[index] building indexes ...", flush=True)
db.execute("CREATE INDEX idx_assoc_study ON gwas_assoc(study_id)")
db.execute("CREATE INDEX idx_assoc_rsid ON gwas_assoc(rsid)")
db.execute("CREATE INDEX idx_studies_cohort ON gwas_studies(cohort, phenotype_category)")
db.commit()
na = db.execute("SELECT count(*) FROM gwas_assoc").fetchone()[0]
ns = db.execute("SELECT count(*) FROM gwas_studies").fetchone()[0]
print(f"[done] {ns:,} studies, {na:,} associations  in {time.time()-t_all:.0f}s; db={os.path.getsize(DB)/1e9:.2f}GB")
db.close()
