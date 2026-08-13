

library(readr)
library(dplyr)

Evo2_score <- read_tsv("Evo2_6Models_Inference_result.6475639.UCSC_Fullanno.LCR_addMAF_Tier.Jan19.txt.gz")

# ================================
# Function to compute correlations (ALL only)
# ================================
run_corr_table_all <- function(df, model_name) {
  
  delta_col  <- paste0(model_name, "_delta_score")
  categories <- c("Q1_0–20%","Q2_20–40%","Q3_40–60%","Q4_60–80%","Q5_80–100%")
  
  stopifnot(delta_col %in% names(df))
  stopifnot("phyloP100way" %in% names(df))
  stopifnot("MAF_tier" %in% names(df))
  stopifnot("Alt_Freq" %in% names(df))
  
  safe_spearman <- function(dat) {
    dat <- dat[is.finite(dat[[delta_col]]) & is.finite(dat$phyloP100way) & is.finite(dat$Alt_Freq), , drop = FALSE]
    if (nrow(dat) <= 2) return(list(n = nrow(dat), rho = NA_real_, p = NA_real_))
    
    sign_factor <- ifelse(dat$Alt_Freq > 0.5, -1, 1)
    ct <- cor.test(dat[[delta_col]] * sign_factor, dat$phyloP100way,
                   method = "spearman", exact = FALSE)
    list(n = nrow(dat), rho = as.numeric(ct$estimate), p = ct$p.value)
  }
  
  out <- data.frame(
    model = character(),
    Category = character(),
    N_variants = integer(),
    rho = numeric(),
    P = numeric(),
    stringsAsFactors = FALSE
  )
  
  # --- all variants ---
  s <- safe_spearman(df)
  out <- rbind(out, data.frame(
    model = model_name,
    Category = "All_variants",
    N_variants = s$n,
    rho = s$rho,
    P = s$p
  ))
  
  # --- per MAF tier ---
  for (cat in categories) {
    sub <- df[df$MAF_tier == cat, , drop = FALSE]
    s2 <- safe_spearman(sub)
    if (s2$n > 2) {
      out <- rbind(out, data.frame(
        model = model_name,
        Category = cat,
        N_variants = s2$n,
        rho = s2$rho,
        P = s2$p
      ))
    }
  }
  
  out[, c("model","Category","N_variants","rho","P")]
}

# Run all models
tbl_7b_noRC      <- run_corr_table_all(Evo2_score, "7b_noRC")
tbl_7b_avgRC     <- run_corr_table_all(Evo2_score, "7b_avgRC")
tbl_7b_weightRC  <- run_corr_table_all(Evo2_score, "7b_weightRC")

tbl_40b_noRC     <- run_corr_table_all(Evo2_score, "40b_noRC")
tbl_40b_avgRC    <- run_corr_table_all(Evo2_score, "40b_avgRC")
tbl_40b_weightRC <- run_corr_table_all(Evo2_score, "40b_weightRC")

all_tables <- rbind(tbl_7b_noRC, tbl_7b_avgRC, tbl_7b_weightRC,
                    tbl_40b_noRC, tbl_40b_avgRC, tbl_40b_weightRC)

write.table(
  all_tables,
  "Evo2_spearman_correlation_summary_phyloP100way_RawScore_MAFbin_ALL.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
