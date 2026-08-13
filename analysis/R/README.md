# Reference R implementations

These are the reference analyses the database reproduces. The Python ports in `analysis/` are
validated against the output tables these scripts produce (see the Validation section of the
top-level README).

The scripts were written against the analysis release of the score table, so their column names
and file paths are those of the source distribution rather than the database schema. They are
included here as the **methodological reference**; run the Python ports for reproducible output.

## `allele_age/` — allele age and recent selection

Uses **site subsampling, not bootstrap resampling**: each of B = 2,000 draws takes 80% of the
sites *without replacement* from every stratum, pools them, and recomputes Spearman ρ. Significance
is read off the 95% subsampling interval (excludes zero → robust), not from a t-test on the ρ
distribution, which would give biased p-values.

| Script | Analysis | Strata | Per-stratum draw |
|---|---|---|---|
| `2a_genomewide_GEVA_subsampling.R` | Genome-wide, Spearman(Δ, GEVA age) on the pruned HapMap3 set, by All / LCR / non-LCR | `AF_bin` | `ceil(0.8n)`, ≥ 1,000, capped at *n* |
| `2d_genomewide_SDS_subsampling.R` | The same, against SDS | `AF_bin` | `ceil(0.8n)`, ≥ 1,000, capped at *n* |
| `2b_HAR_chr_x_AFbin_subsampling.R` | Human Accelerated Regions | `CHR × AF_bin` | `max(1, floor(0.8n))` |
| `2c_RBA_chr_x_AFbin_subsampling.R` | Regulatory-block (ASD) region set | `CHR × AF_bin` | `max(1, floor(0.8n))` |

Note the two per-stratum rules differ between the genome-wide (2a/2d) and region-set (2b/2c)
analyses; the Python ports keep that distinction.

`AF_bin` is equal-width bins of the non-Finnish-European allele frequency
(`AF_0-0.2` … `AF_0.8-1.0`). LCR is defined by the presence of a RepeatMasker `repFamily`
annotation. Allele age is GEVA `AgeMean_Jnt`; selection is `SDS_Final`, re-signed to the
database's allele orientation.

## `conservation/` — conservation and functional stratification

Stratified by `MAF_tier` (NFE minor-allele-frequency quintiles), with bootstrap resampling.

| Script | Analysis |
|---|---|
| `3A_phastCons_by_MAF_tier.R` | Spearman(&#124;Δ&#124;, phastCons100way) overall and per MAF tier |
| `3B_phyloP_signed_delta_by_MAF_tier.R` | Spearman(Δ × minor-allele sign, phyloP100way) — the **raw signed** Δ, not its absolute value |
| `3G_cCRE_category.R` | &#124;Δ&#124; by ENCODE cCRE category, Mann-Whitney against the *Intergenic* reference |
| `3H_coding_consequence.R` | &#124;Δ&#124; by coding consequence, Mann-Whitney against all other consequences |
