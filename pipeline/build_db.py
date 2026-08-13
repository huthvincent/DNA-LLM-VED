#!/usr/bin/env python3
"""
Build a Datasette-ready SQLite database from the re-annotated Evo2 inference TSV.

Source : Jun18/Evo2_6Models_Inference_result.6475578.UCSC_Fullanno.LCR_gnomAD.joint.Frequency.June18.txt.gz
Output : evo2.db   (table: variants  + FTS + indexes)

The June-18 re-annotation (re_annotate.R) merged gnomAD v4.1 *joint* allele
frequencies onto the variant set, renamed every column to a human-readable
scheme, dropped the old single-cohort frequency columns (Alt_Freq / Callrate /
MAF / MAF_tier) and added joint + 8 ancestry-specific AFs. The file already
ships final column names and a ready Variant_ID (= Chromosome:Position:Ref:Alt),
so the build just streams it in with type inference, NULL handling, indexes and
an FTS5 index over RSID / Gene / amino-acid change.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
import config

import csv, gzip, sqlite3, sys, time

SRC = str(config.VARIANT_TSV)
DB  = str(config.DB_PATH)

NULL_TOKENS = {"", ".", "NA", "na", "N/A", "NaN", "nan", "None", "null"}
BATCH = 20000
csv.field_size_limit(1 << 24)  # some annotation columns can be long


def clean(col: str) -> str:
    """The June-18 header names are already valid identifiers; just normalise."""
    return col.strip().replace(".", "_").replace("-", "_").replace(" ", "_")


def col_type(name: str) -> str:
    """INTEGER for Position; REAL for allele frequencies, Evo2 scores and
    conservation; TEXT for everything else."""
    if name == "Position":
        return "INTEGER"
    if name.endswith("_AF") or name.startswith("Evo2_") or name in ("PhastCons_100way", "PhyloP_100way"):
        return "REAL"
    return "TEXT"


def main():
    t0 = time.time()
    with gzip.open(SRC, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        raw_header = next(reader)
        cols = [clean(c) for c in raw_header]

        # de-duplicate any repeated cleaned names defensively
        seen, final_cols = {}, []
        for c in cols:
            if c in seen:
                seen[c] += 1; final_cols.append(f"{c}_{seen[c]}")
            else:
                seen[c] = 0; final_cols.append(c)
        cols = final_cols
        types = {c: col_type(c) for c in cols}

        real_idx = {i for i, c in enumerate(cols) if types[c] == "REAL"}
        int_idx  = {i for i, c in enumerate(cols) if types[c] == "INTEGER"}
        ncol = len(cols)

        print(f"[schema] {ncol} columns")
        for c in cols:
            print(f"    {c:28s} {types[c]}")
        sys.stdout.flush()

        db = sqlite3.connect(DB)
        db.executescript(
            "PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;"
            "PRAGMA temp_store=MEMORY; PRAGMA cache_size=-400000;"
        )
        coldefs = ",\n  ".join(f'"{c}" {types[c]}' for c in cols)
        db.execute("DROP TABLE IF EXISTS variants")
        db.execute(f'CREATE TABLE variants (\n  {coldefs}\n)')

        placeholders = ",".join(["?"] * ncol)
        insert = f"INSERT INTO variants VALUES ({placeholders})"

        def conv(row):
            if len(row) != ncol:
                row = (row + [None] * ncol)[:ncol]
            out = [None] * ncol
            for i in range(ncol):
                v = row[i]
                if v is None or v in NULL_TOKENS:
                    out[i] = None
                elif i in real_idx:
                    try: out[i] = float(v)
                    except ValueError: out[i] = None
                elif i in int_idx:
                    try: out[i] = int(float(v))
                    except ValueError: out[i] = None
                else:
                    out[i] = v
            return out

        buf, n, truncated = [], 0, False
        try:
            for row in reader:
                buf.append(conv(row))
                if len(buf) >= BATCH:
                    db.executemany(insert, buf); buf.clear()
                    n += BATCH
                    if n % 500000 == 0:
                        print(f"[load] {n:,} rows  {time.time()-t0:,.0f}s"); sys.stdout.flush()
        except (EOFError, OSError) as e:
            truncated = True
            print(f"[warn] source gzip truncated near row {n+len(buf):,}: {e}")
        if buf:
            db.executemany(insert, buf); n += len(buf)
        db.commit()
        print(f"[load] DONE {n:,} rows  truncated={truncated}  {time.time()-t0:,.0f}s"); sys.stdout.flush()

    # indexes on the columns users actually query / facet / sort by
    idx_cols = [
        ("idx_chr_pos", '"Chromosome","Position"'),
        ("idx_rsid", '"RSID"'),
        ("idx_gene", '"Gene"'),
        ("idx_variant", '"Variant_ID"'),
        ("idx_7b_delta", '"Evo2_7B_NoRC_Delta"'),
        ("idx_40b_delta", '"Evo2_40B_NoRC_Delta"'),
        ("idx_func", '"Functional_Annotation"'),
        ("idx_exonic", '"Exonic_Function"'),
        ("idx_ccre", '"ENCODE_cCRE_Annotation"'),
        ("idx_repclass", '"Repeat_Class"'),
        ("idx_nfe_af", '"Non_Finnish_European_AF"'),
        ("idx_joint_af", '"Joint_AF"'),
    ]
    for name, expr in idx_cols:
        try:
            db.execute(f"CREATE INDEX {name} ON variants({expr})")
            print(f"[index] {name}"); sys.stdout.flush()
        except sqlite3.OperationalError as e:
            print(f"[index] SKIP {name}: {e}")

    # FTS5 free-text search over identifiers / gene / protein change
    fts_cols = [c for c in ("RSID", "Gene", "Amino_Acid_Change") if c in cols]
    if fts_cols:
        cl = ", ".join(f'"{c}"' for c in fts_cols)
        db.execute(
            f"CREATE VIRTUAL TABLE variants_fts USING fts5({cl}, "
            f"content='variants', content_rowid='rowid')"
        )
        db.execute(f"INSERT INTO variants_fts(rowid, {cl}) SELECT rowid, {cl} FROM variants")
        print(f"[fts] indexed {fts_cols}"); sys.stdout.flush()

    print("[analyze] ..."); sys.stdout.flush()
    db.execute("ANALYZE"); db.commit()
    db.execute("VACUUM"); db.close()
    print(f"[done] {time.time()-t0:,.0f}s total")


if __name__ == "__main__":
    main()
