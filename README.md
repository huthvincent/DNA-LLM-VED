# Evo2VED

**A database of DNA large language model–derived scores for genome-wide variant interpretation.**

Evo2VED provides zero-shot delta-likelihood scores from the [Evo 2](https://github.com/ArcInstitute/evo2)
DNA foundation model for **6,475,578 common human variants** (GRCh38, autosomes, gnomAD
non-Finnish-European MAF ≥ 5%, predominantly non-coding), integrated with Ensembl, RepeatMasker,
ENCODE cCRE, phastCons/phyloP, gnomAD v4.1 joint allele frequencies, allele age (GEVA), recent
selection (SDS), and GWAS/S-LDSC results.

**Live database:** <https://huthvincent-evo2-database.hf.space>

This repository contains the **code** — the database build pipeline, the statistical analyses, the
web application, and the figure scripts. It does not contain the variant data itself (see
[Data](#data)).

---

## Scores

Each variant carries a delta score `Δ = log P(ALT) − log P(REF)` under **six configurations**:
two model scales × three reverse-complement (RC) strategies.

| | no RC | average RC | weighted RC |
|---|---|---|---|
| **Evo2-7B** | `Evo2_7B_NoRC_Delta` | `Evo2_7B_AvgRC_Delta` | `Evo2_7B_WeightedRC_Delta` |
| **Evo2-40B** | `Evo2_40B_NoRC_Delta` | `Evo2_40B_AvgRC_Delta` | `Evo2_40B_WeightedRC_Delta` |

`noRC` scores the forward strand only; `avgRC` averages forward and reverse-complement scores;
`weightedRC` weights each strand by its per-strand likelihood. The underlying REF and ALT
log-likelihoods are stored alongside every delta.

## Analysis modules

The web application exposes four interactive modules, all computed server-side:

| Module | Route | What it does |
|---|---|---|
| Regional analysis | `/viewer` | Δ stratified by functional annotation and coding consequence within a gene or interval; 6-configuration concordance; locus landscape |
| Region-set analysis | `/regionset` | Conservation, allele age, allele-frequency and functional stratification for a built-in region set (genome-wide, HAR) or an uploaded BED |
| GWAS / heritability | `/gwas` | S-LDSC coefficient z-scores across UK Biobank, FinnGen and MVP traits, ranked within and across phenotype categories |
| Co-localization | `/coloc` | Genome-browser view lining up the Evo2 Δ track, GWAS −log₁₀(P) tracks and variant annotation over a region |

## Statistical methods

The analyses reproduce the reference R implementations in [`analysis/R/`](analysis/R); the Python
ports are validated against the published reference outputs (see [Validation](#validation)).

**Conservation (Tasks 4A / 4B)** — stratified by `MAF_tier` (NFE minor-allele-frequency quintiles),
with bootstrap resampling.

- *4A phastCons*: Spearman(|Δ|, phastCons100way)
- *4B phyloP*: Spearman(Δ × minor-allele sign, phyloP100way) — the **raw signed** Δ, where the sign
  is −1 when the NFE allele frequency exceeds 0.5.

**Allele age and selection (Tasks 5 / 6)** — **site subsampling, not bootstrap resampling.** In each
of *B* = 2,000 draws, 80% of the sites are taken **without replacement** from every stratum
(`AF_bin` genome-wide, `CHR × AF_bin` for region sets), pooled, and Spearman ρ is recomputed.

- The reported **ρ is the full-data point estimate**.
- The interval is the 2.5/97.5 percentile of the subsampling distribution.
- **Significance is judged by whether that interval excludes zero.**
- The directional rate `mean(ρ ≤ 0)` is floored at 1/(B+1) and is reported for reference only —
  it is *not* a conventional p-value.

Metrics: allele age is GEVA `AgeMean_Jnt`; selection is `SDS_Final`, the Singleton Density Score
re-signed to match this database's allele orientation. Low-complexity regions (LCR) are defined by
the presence of a RepeatMasker `repFamily` annotation.

**Functional stratification (Tasks 7G / 7H)** — |Δ| by ENCODE cCRE category (Mann-Whitney against
the *Intergenic* reference) and by coding consequence (against all other consequences), both
Benjamini-Hochberg adjusted.

## Repository layout

```
config.py            central path configuration — every script resolves paths through this
pipeline/            database construction and ingest
  build_db.py                 build evo2.db from the annotated score table
  add_maf_tier.py             derive MAF_tier (NFE minor-allele quintiles)
  ingest_allele_age_v2.py     allele_age table (GEVA AgeMean_Jnt + SDS_Final) and AF_bin
  ingest_module2.py           region_sets (HAR) and the legacy allele-age ingest
  ingest_gwas_p_v2.py         GWAS associations from the pre-filtered P ≤ 1e-4 release
  ingest_gwas_p.py            legacy ingest for the streamed-ZIP summary statistics
  build_gwas.py               normalise the S-LDSC result files into gwas_ldsc.json
  fix_gwas_annotation.py      re-map trait names/categories from the authoritative annotation
  verify_db.py                post-build sanity checks
analysis/            statistics
  precompute_regionset.py     genome-wide results for all 6 models -> regionset_genomewide.json
  precompute_stats.py         genome-wide reference means -> genome_stats.json
  validate_subsampling.py     validate the genome-wide subsampling against the reference
  validate_har.py             validate the CHR x AF_bin region-set subsampling against the reference
  R/                          reference R implementations (conservation, allele age)
webapp/              Datasette application
  plugins/analysis.py         server-side compute routes (the /api/... endpoints)
  templates/, static/         front end
  metadata.yml                Datasette configuration
  Dockerfile, entrypoint.sh   container image
figures/             manuscript figure scripts
examples/            example inputs (HAR region BED)
docs/                deployment notes and supplementary material
```

## Data

The variant table is **not** in this repository — it is several gigabytes. The pipeline expects the
source files under a data root you control:

```bash
export EVO2VED_DATA=/path/to/evo2ved-data
python config.py            # print every resolved path
```

`config.py` documents each expected input and the environment variables that override it
(`EVO2VED_DB`, `EVO2VED_SOURCE`, `EVO2VED_OUT`, `EVO2VED_FIGURES`). Data access is described on the
live site; please contact the corresponding authors for the source tables.

## Running

```bash
pip install -r requirements.txt

# 1. build the database (long-running)
python pipeline/build_db.py
python pipeline/add_maf_tier.py
python pipeline/ingest_allele_age_v2.py
python pipeline/ingest_module2.py
python pipeline/ingest_gwas_p_v2.py
python pipeline/build_gwas.py
python pipeline/verify_db.py

# 2. pre-compute the genome-wide statistics served by /regionset
python analysis/precompute_stats.py
python analysis/precompute_regionset.py

# 3. serve
cd webapp && ./entrypoint.sh          # or: docker build -t evo2ved . && docker run -p 7860:7860 evo2ved
```

The genome-wide region-set results are pre-computed offline for all six models so the
`/regionset` page can display and export any model instantly, with no per-request computation.

## Validation

`analysis/validate_subsampling.py` and `analysis/validate_har.py` re-run the subsampling procedure
against the reference outputs distributed with the R code. Both reproduce the reference exactly
where the procedure is deterministic:

| Check | Result |
|---|---|
| Genome-wide sites per draw (All / LCR / non-LCR) | 40,863 / 13,244 / 27,806 — exact |
| Genome-wide ρ and interval | within 2×10⁻⁴ (Monte-Carlo noise only) |
| HAR sites per stratum and per draw | exact (1,519 total; 1,180 per draw) |
| HAR full-data ρ | matches to 1×10⁻¹⁷ (machine precision) |

Remaining differences are the RNG only: R's `sample()` and NumPy's generator draw different
subsets, so percentile bounds differ in the fourth decimal.

## Contributors

- **Rui Zhu** ([@huthvincent](https://github.com/huthvincent)) — Yale School of Medicine
- **Xiaopu Zhou** — The Hospital for Sick Children (SickKids)

## Citation

Manuscript in preparation. Please cite the live database until the paper appears:

> Evo2VED: a database of DNA large language model–derived scores for genome-wide variant
> interpretation. <https://huthvincent-evo2-database.hf.space>

Supplementary material for the manuscript is in [`docs/supplementary/`](docs/supplementary).

## License

Released under the [MIT License](LICENSE). The variant data is distributed separately and is not
covered by this license.
