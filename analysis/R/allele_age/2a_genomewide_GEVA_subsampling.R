#setwd("C:\\Users\\fredz\\OneDrive\\2026Jan\\EVO2_noncoding_Fix\\Analysis\\Figure2\\2a")


library(readr)
library(dplyr)
library(data.table)

setDTthreads(0)

# -----------------------------
# 1) Load data
# -----------------------------
Evo2_score <- read_tsv(
  "Evo2_6Models_Inference_result.GEVA_Fullanno.LCR_addAF_Tier_4629351.July12.txt.gz",
  show_col_types = FALSE
) %>% as.data.table()

Evo2_score[!is.na(Repeat_Family), Is_LCR := "LCR"]
Evo2_score[ is.na(Repeat_Family), Is_LCR := "nonLCR"]

# pruned HM3 list
prune <- fread(
  "1000GP3_EUR_MAF5_gAD-call0.8_968159_hm3_allele_age.pruned.Final.prune.in",
  header = FALSE
)
setnames(prune, "V1", "rsid_pruned")

Evo2_score_prune <- Evo2_score[RSID %in% prune$rsid_pruned]

# -----------------------------
# 2) Bootstrap correlation (stratified by af_tier)
# -----------------------------
af_levels <-  c("AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0")
B <- 2000
n_total <- 5000
n_per_maf <- n_total / length(af_levels)  # = 1000

stopifnot(n_per_maf == as.integer(n_per_maf))

bootstrap_spearman_fixedFrac <- function(df, model_name,
                                         region_col = "Is_LCR",
                                         region_levels = c("All","LCR","nonLCR"),
                                         af_col = "AF_bin",
                                         af_levels = NULL,
                                         B = 500,
                                         frac_by_region = c(All=0.10, LCR=0.10, nonLCR=0.10),
                                         seed = 1,
                                         with_replacement = FALSE,
                                         rounding = c("ceiling","floor","round"),
                                         min_per_bin = 1) {
  
  if (is.null(af_levels)) {
    af_levels <- c("AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0")
  }
  rounding <- match.arg(rounding)
  
  delta_col <- paste0(model_name, "_Delta")
  stopifnot(delta_col %in% names(df))
  stopifnot("AgeMedian_Jnt" %in% names(df))
  stopifnot(af_col %in% names(df))
  stopifnot(all(region_levels %in% names(frac_by_region)))
  
  set.seed(seed)
  
  # Pre-filter
  df2 <- df[
    is.finite(get(delta_col)) & is.finite(AgeMedian_Jnt) & !is.na(get(af_col))
  ]
  
  out <- data.table(
    model = character(),
    Region = character(),
    bootstrap_id = integer(),
    N_variants = integer(),
    rho = numeric(),
    P = numeric()
  )
  
  calc_spearman <- function(x, y) {
    rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
    p   <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
    list(rho = as.numeric(rho), p = as.numeric(p))
  }
  
  # convert fraction -> integer sample size for one bin
  frac_to_n <- function(frac, n_avail) {
    raw <- frac * n_avail
    n <- switch(rounding,
                ceiling = ceiling(raw),
                floor   = floor(raw),
                round   = as.integer(round(raw)))
    n <- max(n, min_per_bin)
    if (!with_replacement) n <- min(n, n_avail)
    n
  }
  
  for (reg in region_levels) {
    dreg <- if (reg == "All") df2 else df2[get(region_col) == reg]
    if (nrow(dreg) == 0) next
    
    frac <- as.numeric(frac_by_region[[reg]])
    
    draw_once <- function() {
      idx_all <- integer(0)
      
      for (lvl in af_levels) {
        sub <- which(dreg[[af_col]] == lvl)
        n_avail <- length(sub)
        if (n_avail == 0) next  # bin absent -> contributes 0
        
        k <- frac_to_n(frac, n_avail)
        if (k <= 0) next
        
        idx <- sample(sub, size = k, replace = with_replacement)
        idx_all <- c(idx_all, idx)
      }
      idx_all
    }
    
    for (b in seq_len(B)) {
      idx <- draw_once()
      if (length(idx) <= 2) next
      
      x <- dreg[[delta_col]][idx]
      y <- dreg[["AgeMedian_Jnt"]][idx]
      sp <- calc_spearman(x, y)
      
      out <- rbind(out, data.table(
        model = model_name,
        Region = reg,
        bootstrap_id = b,
        N_variants = length(idx),
        rho = sp$rho,
        P = sp$p
      ))
    }
  }
  
  out[]
}


# Run bootstrap for one model
af_levels <- c("AF_0-0.2", "AF_0.2-0.4", "AF_0.4-0.6", "AF_0.6-0.8", "AF_0.8-1.0")

tbl_boot <- bootstrap_spearman_fixedFrac(
  df = Evo2_score_prune,
  model_name = "Evo2_40B_AvgRC",
  B = 2000,
  frac_by_region = c(All=0.8, LCR=0.8, nonLCR=0.8),
  seed = 123,
  with_replacement = FALSE,     # recommended for bootstrap
  rounding = "ceiling",
  min_per_bin = 1000
)

# -----------------------------
# 3) Summarize bootstrap CI per region
# -----------------------------
ci_tbl <- tbl_boot[, {
  r <- rho[is.finite(rho)]
  Bn <- length(r)
  
  list(
    B = Bn,
    N_each = unique(N_variants),
    rho_mean = mean(r),
    rho_median = median(r),
    
    # variability
    rho_sd  = if (Bn > 1) sd(r) else NA_real_,
    rho_sem = if (Bn > 1) sd(r) / sqrt(Bn) else NA_real_,
    
    # 95% bootstrap CI
    rho_ci2.5  = quantile(r, 0.025),
    rho_ci97.5 = quantile(r, 0.975),
    
    # one-tailed bootstrap test: H1 = rho > 0
    prop_gt0 = mean(r > 0),
    p_one_tailed_gt0 = mean(r <= 0)
  )
}, by = .(model, Region)]



# -----------------------------
# 4) Write outputs
# -----------------------------
fwrite(tbl_boot,
       file = "Bootstrap2000_Spearman_Evo2_40B_AvgRC_Delta_All_LCR_nonLCR_strataf_80Percent_sites_AF_bin.tsv",
       sep = "\t")

fwrite(ci_tbl,
       file = "Bootstrap2000_Spearman_Evo2_40B_AvgRC_Delta_All_LCR_nonLCR_strataf_80Percent_sites_AF_bin.CI.tsv",
       sep = "\t")
