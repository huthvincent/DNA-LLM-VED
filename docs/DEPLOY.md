# Evo2 Variant Effect Database — deployment runbook

Everything below assumes you are in `/Users/rui/Desktop/papers/NAR/database/evo2_database`
and have already built `evo2.db` (3.4 GB) and `evo2.parquet` (796 MB) with `build_db.py`.

```
evo2_database/
├── evo2.db              # SQLite (Datasette serves this)
├── evo2.parquet         # columnar, for download / DuckDB / HF viewer
├── metadata.yml         # Datasette config (titles, facets, canned queries)
├── build_db.py          # rebuild the DB from the source .gz
├── verify_db.py         # sanity checks
├── space/               # HuggingFace Docker Space (Datasette)
│   ├── Dockerfile  ├── entrypoint.sh  ├── metadata.yml  └── README.md
└── dataset/README.md    # HuggingFace dataset card + data dictionary
```

> **Before the public release:** the source `.gz` was truncated (last ~60–166 of
> 6,475,639 rows lost). Re-export a clean `.gz`, re-run `python3 build_db.py`, and
> re-export Parquet to ship the complete set. The current build (6,475,578 rows) is
> fine for the demo and the pre-submission inquiry.

---

## Path A — HuggingFace (recommended: free, no credit card, you already have an account)

**Backend = a HF *dataset* repo holding the data; front-end = a HF *Space* running Datasette.**

```bash
pip install -U "huggingface_hub[cli]"
huggingface-cli login          # paste a WRITE token from huggingface.co/settings/tokens

DATA=huthvincent/Evo2-Variant-DB      # data repo (or reuse Evo2_7b_scores)
SPACE=huthvincent/evo2-database       # the web app

# 1. data repo: upload Parquet (+ optional SQLite) + card
huggingface-cli repo create "$DATA" --type dataset -y
huggingface-cli upload "$DATA" evo2.parquet      evo2.parquet --repo-type dataset
huggingface-cli upload "$DATA" evo2.db           evo2.db      --repo-type dataset
huggingface-cli upload "$DATA" dataset/README.md README.md    --repo-type dataset

# 2. Datasette Space (Docker SDK)
huggingface-cli repo create "$SPACE" --type space --space_sdk docker -y
huggingface-cli upload "$SPACE" space/ . --repo-type space

# 3. tell the Space where to fetch the DB on startup
python3 - <<'PY'
from huggingface_hub import add_space_variable
add_space_variable("huthvincent/evo2-database", "DB_URL",
  "https://huggingface.co/datasets/huthvincent/Evo2-Variant-DB/resolve/main/evo2.db")
PY
```

Live at `https://huggingface.co/spaces/huthvincent/evo2-database` (Datasette UI + JSON API).
Uploading the Parquet also lights up HF's own **Dataset Viewer + SQL console** for free.

> **Cold-start note (free CPU Space):** ephemeral disk ⇒ the 3.4 GB DB re-downloads on
> each restart. To remove that wait, either enable **persistent storage** in Space settings,
> or serve `evo2.parquet` via the `datasette-parquet` plugin (smaller, DuckDB-backed).

---

## Path B — Google Cloud Run (one command, scales to zero ≈ $0 idle, fastest cold start)

Needs `gcloud` installed + a GCP project (free tier covers a low-traffic DB).
Datasette has this built in — the DB is baked into the image, so **no re-download**:

```bash
pip install datasette
gcloud auth login && gcloud config set project YOUR_PROJECT
datasette publish cloudrun evo2.db \
  -m metadata.yml \
  --service evo2-db \
  --memory 4Gi \
  --extra-options "--setting sql_time_limit_ms 12000 --setting max_returned_rows 1000"
```

Gives a stable `https://evo2-db-xxxxx.run.app` URL.

---

## Make it NAR-grade

- **Custom domain** (~$10/yr): map e.g. `evo2db.org` to the Space/Cloud Run service so the
  URL is permanent regardless of host. NAR reviewers weigh URL persistence.
- **Zenodo DOI** (free): deposit `evo2.parquet` for a citable, archived snapshot; link the DOI
  from the dataset card and the paper.
- **No login wall** — already satisfied (public, anonymous read). This is a hard NAR requirement.

## Local preview

```bash
./.venv/bin/datasette -i evo2.db -m metadata.yml --port 8765
# → http://127.0.0.1:8765
```
