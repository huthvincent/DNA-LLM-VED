# syntax=docker/dockerfile:1
# Datasette server for the Evo2 Variant Effect Database, packaged for a
# HuggingFace Docker Space (or any container host). Listens on $PORT (7860).
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir "datasette==0.65.2"

WORKDIR /app

# Bake the database into the image FIRST, so this heavy layer stays cached and
# front-end changes (templates / metadata / entrypoint) rebuild in seconds
# without re-downloading the 3.4 GB db.
# The dataset repo is PRIVATE, so authenticate with the HF_TOKEN build secret if
# present (set under Space Settings -> Variables and secrets -> HF_TOKEN).
# Falls back to an unauthenticated fetch if no secret is set (e.g. public dataset).
ARG DB_URL="https://huggingface.co/datasets/huthvincent/Evo2-Variant-DB/resolve/main/evo2.db"
# DBVER busts the build cache for the fetch layer whenever the dataset's evo2.db changes.
# Bump this every time you re-upload a new evo2.db, or the old db stays cached/baked.
ARG DBVER="2026-07-24-alleleage-v2-afbin-gwasanno"
RUN --mount=type=secret,id=HF_TOKEN,mode=0444,required=false \
    echo "fetching db version: $DBVER"; \
    T=""; if [ -f /run/secrets/HF_TOKEN ]; then T="$(cat /run/secrets/HF_TOKEN)"; fi; \
    if [ -n "$T" ]; then \
      wget -q --header="Authorization: Bearer $T" -O /app/evo2.db "$DB_URL"; \
    else \
      wget -q -O /app/evo2.db "$DB_URL"; \
    fi; \
    echo "baked db [$DBVER]: $(du -h /app/evo2.db | cut -f1)"

# Compute dependencies for the analysis plugin (numpy / scipy / pandas).
# Placed AFTER the heavy db layer so changing it never re-downloads the 3.4 GB db.
RUN pip install --no-cache-dir numpy scipy pandas

# Front-end + config + compute plugin (cheap, changes often -> after the heavy layers)
COPY metadata.yml entrypoint.sh genome_stats.json gwas_ldsc.json regionset_genomewide.json ./
COPY templates/ ./templates/
COPY static/ ./static/
COPY plugins/ ./plugins/
RUN chmod +x entrypoint.sh

ENV PORT=7860
EXPOSE 7860
CMD ["./entrypoint.sh"]
