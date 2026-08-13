"""Central path configuration for the Evo2VED pipeline.

Every script in this repository resolves its inputs and outputs through this module,
so nothing is hard-coded to one machine. Override any location with an environment
variable; otherwise everything lives under `data/` next to this file.

    EVO2VED_DATA     root for all data (default: <repo>/data)
    EVO2VED_DB       SQLite database        (default: $EVO2VED_DATA/evo2.db)
    EVO2VED_SOURCE   raw source annotations (default: $EVO2VED_DATA/source)
    EVO2VED_OUT      generated artefacts    (default: <repo>/webapp, so the app picks them up)
    EVO2VED_FIGURES  figure output          (default: <repo>/figures/output)

Typical use:

    export EVO2VED_DATA=/scratch/evo2ved
    python pipeline/build_db.py

Nothing here creates or downloads data; see README.md for how to obtain the inputs.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent


def _p(env: str, default: Path) -> Path:
    return Path(os.environ.get(env, default)).expanduser()


# --- roots -----------------------------------------------------------------
DATA_ROOT = _p("EVO2VED_DATA", REPO_ROOT / "data")
SOURCE_ROOT = _p("EVO2VED_SOURCE", DATA_ROOT / "source")
OUT_ROOT = _p("EVO2VED_OUT", REPO_ROOT / "webapp")
FIGURE_DIR = _p("EVO2VED_FIGURES", REPO_ROOT / "figures" / "output")

# --- database --------------------------------------------------------------
DB_PATH = _p("EVO2VED_DB", DATA_ROOT / "evo2.db")

# --- raw inputs (see README for provenance) --------------------------------
# Evo2 scores + UCSC/gnomAD annotation, one row per variant
VARIANT_TSV = SOURCE_ROOT / "Evo2_6Models_Inference_result.UCSC_Fullanno.txt.gz"
# Allele-age / selection tables (July-12 harmonised release)
GEVA_TSV = SOURCE_ROOT / "allele_age" / "Evo2_6Models_Inference_result.GEVA.txt.gz"
SDS_TSV = SOURCE_ROOT / "allele_age" / "Evo2_6Models_Inference_result.SDS.txt.gz"
PRUNE_IN = SOURCE_ROOT / "allele_age" / "hm3_allele_age.pruned.Final.prune.in"
# Region sets (BED-like, chrom/start/end)
HAR_BED = SOURCE_ROOT / "regions" / "HAR_hg38_bed_region.txt"
# S-LDSC results and the authoritative per-cohort trait annotation
LDSC_DIR = SOURCE_ROOT / "ldsc" / "LDSC_results"
GWAS_P_DIR = SOURCE_ROOT / "gwas" / "LDSC_P_1E-4"
GWAS_ANNO_DIR = SOURCE_ROOT / "gwas_annotation"

# --- generated artefacts (baked into the web app) --------------------------
GENOME_STATS_JSON = OUT_ROOT / "genome_stats.json"
GWAS_LDSC_JSON = OUT_ROOT / "gwas_ldsc.json"
REGIONSET_GW_JSON = OUT_ROOT / "regionset_genomewide.json"
MAF_TIER_TSV = DATA_ROOT / "MAF_tier_annotation.tsv.gz"
CASE_DATA_JSON = FIGURE_DIR / "case_data.json"
SCREENS_DIR = FIGURE_DIR / "screens"

# --- local server used by the figure/screenshot scripts --------------------
LOCAL_SERVER = os.environ.get("EVO2VED_SERVER", "http://127.0.0.1:8765")


def ensure_dirs() -> None:
    """Create the output directories a script is about to write into."""
    for d in (DATA_ROOT, OUT_ROOT, FIGURE_DIR):
        d.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    for k, v in sorted(globals().items()):
        if k.isupper() and isinstance(v, (Path, str)):
            print(f"{k:20s} {v}")
