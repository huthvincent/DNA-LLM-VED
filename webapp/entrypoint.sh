#!/usr/bin/env bash
# Serve evo2.db with Datasette. If the DB is not already present (e.g. mounted
# for local testing), download it once from $DB_URL (a HuggingFace resolve URL).
set -euo pipefail

DB="${DB_PATH:-/app/evo2.db}"

if [ ! -f "$DB" ]; then
  : "${DB_URL:?Set DB_URL to the evo2.db resolve URL (or mount the DB at /app/evo2.db)}"
  echo "[entrypoint] downloading database from: $DB_URL"
  wget -q -O "$DB" "$DB_URL"
  echo "[entrypoint] downloaded $(du -h "$DB" | cut -f1)"
fi

echo "[entrypoint] serving $DB on port ${PORT:-7860}"
exec datasette -i "$DB" -m metadata.yml \
  --template-dir templates --static static:static \
  --plugins-dir plugins \
  --host 0.0.0.0 --port "${PORT:-7860}" \
  --setting sql_time_limit_ms 20000 \
  --setting facet_time_limit_ms 6000 \
  --setting max_returned_rows 6000 \
  --setting default_facet_size 20 \
  --setting suggest_facets off \
  --setting allow_download off \
  --setting allow_csv_stream off \
  --setting max_csv_mb 5 \
  --cors
