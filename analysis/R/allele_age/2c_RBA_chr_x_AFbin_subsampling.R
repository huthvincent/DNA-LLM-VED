setwd("G:\\Project\\Evo2VEP\\Allele_age\\2c")

library(data.table)
library(readr)
library(dplyr)

## =========================
## 0) Read data
## =========================
Evo2_score <- as.data.table(
  read_tsv("../2a/Evo2_6Models_Inference_result.Alelle_Age.UCSC_Fullanno.LCR_addAF_Tier_4629351.July12.txt.gz",
           show_col_types = FALSE)
)

Evo2_score[, Is_LCR := fifelse(!is.na(Repeat_Family), "LCR", "nonLCR")]

RBA <- fread("RBA_ASD.txt", header = FALSE)
RBA=RBA[-1,]
setnames(RBA, c("CHR", "start", "end"))
regions_dt <- as.data.table(RBA)

## normalize CHR format
Evo2_score[, CHR := as.character(Chromosome)]
Evo2_score[!grepl("^chr", CHR), CHR := paste0("chr", CHR)]
regions_dt[, CHR := as.character(CHR)]
regions_dt[!grepl("^chr", CHR), CHR := paste0("chr", CHR)]

## =========================
## 1) Region overlap (foverlaps)
## =========================
Evo2_score[, vid := .I]
vars_iv <- Evo2_score[, .(CHR, start = as.integer(Position), end = as.integer(Position), vid)]

regs_iv <- regions_dt[, .(CHR, start = as.integer(start), end = as.integer(end))]
regs_iv[, region_id := paste0(CHR, ":", start, "-", end)]

setkey(regs_iv, CHR, start, end)
setkey(vars_iv, CHR, start, end)

hits <- foverlaps(vars_iv, regs_iv, type = "within", nomatch = 0L)

anno <- hits[, .(
  in_region   = TRUE,
  n_regions   = .N,
  regions_hit = paste(unique(region_id), collapse = ";")
), by = vid]

variants_annotated <- merge(Evo2_score, anno, by = "vid", all.x = TRUE)

variants_annotated[is.na(in_region), `:=`(
  in_region = FALSE,
  n_regions = 0L,
  regions_hit = NA_character_
)]
variants_annotated[, vid := NULL]

## keep only in-region variants (LATEST choice)
variants_annotated_clean <- variants_annotated[in_region == TRUE]

## quick plot
plot(variants_annotated_clean$AgeMean_Jnt,
     variants_annotated_clean$`Evo2_40B_AvgRC_Delta`)


## ============================================================
## 3) Stratified site subsampling
##
## In each iteration:
##   - retain every chromosome;
##   - within each CHR × AF_bin stratum;
##   - sample 80% of sites without replacement;
##   - calculate Spearman rho for all sites and each AF bin.
## ============================================================

boot_by_chr_maf_site_spearman <- function(
    df,
    B = 2000,
    seed = 1,
    alternative = c("greater", "less"),
    prop_site = 0.8,
    return_folds = TRUE
) {
  
  alternative <- match.arg(alternative)
  
  if (prop_site <= 0 || prop_site > 1) {
    stop("prop_site must be greater than 0 and no greater than 1.")
  }
  
  set.seed(seed)
  
  ## ----------------------------------------------------------
  ## Prepare analysis data
  ## ----------------------------------------------------------
  
  d <- df %>%
    dplyr::filter(
      !is.na(CHR),
      !is.na(AF_bin)
    ) %>%
    dplyr::transmute(
      CHR = as.character(CHR),
      AF_bin = as.character(AF_bin),
      AgeMean_Jnt = as.numeric(AgeMean_Jnt),
      score = as.numeric(`Evo2_40B_AvgRC_Delta`)
    ) %>%
    dplyr::filter(
      stats::complete.cases(.)
    )
  
  if (nrow(d) < 3L) {
    stop("Too few complete observations.")
  }
  
  ## Convert to data.table for efficient stratified sampling
  d <- data.table::as.data.table(d)
  d[, row_id := .I]
  
  all_chrs <- sort(unique(d$CHR))
  n_chr_total <- length(all_chrs)
  
  if (n_chr_total < 1L) {
    stop("No chromosomes with complete observations were found.")
  }
  
  ## ----------------------------------------------------------
  ## Spearman correlation helper
  ## ----------------------------------------------------------
  
  point_rho <- function(sub) {
    
    if (nrow(sub) < 3L) {
      return(NA_real_)
    }
    
    if (
      data.table::uniqueN(sub$AgeMean_Jnt) < 2L ||
      data.table::uniqueN(sub$score) < 2L
    ) {
      return(NA_real_)
    }
    
    suppressWarnings(
      stats::cor(
        sub$AgeMean_Jnt,
        sub$score,
        method = "spearman",
        use = "complete.obs"
      )
    )
  }
  
  ## ----------------------------------------------------------
  ## Full-data point estimates
  ## ----------------------------------------------------------
  
  rho_all_hat <- point_rho(d)
  
  rho_bin_hat <- d[
    ,
    .(
      rho = point_rho(.SD),
      N = .N,
      n_chr = uniqueN(CHR)
    ),
    by = AF_bin
  ]
  
  ## ----------------------------------------------------------
  ## Define CHR × AF-bin strata
  ## ----------------------------------------------------------
  
  stratum_key <- interaction(
    d$CHR,
    d$AF_bin,
    drop = TRUE,
    lex.order = TRUE
  )
  
  stratum_indices <- split(
    d$row_id,
    stratum_key
  )
  
  ## Number sampled from each CHR × AF-bin stratum
  ##
  ## floor() ensures no more than 80% are selected.
  ## max(1, ...) ensures that small non-empty strata remain
  ## represented in every iteration.
  stratum_sample_sizes <- vapply(
    stratum_indices,
    function(index_vector) {
      
      n_available <- length(index_vector)
      
      max(
        1L,
        as.integer(
          floor(prop_site * n_available)
        )
      )
    },
    integer(1)
  )
  
  ## ----------------------------------------------------------
  ## Storage
  ## ----------------------------------------------------------
  
  rhos_all <- rep(
    NA_real_,
    B
  )
  
  bins <- sort(
    unique(d$AF_bin)
  )
  
  rhos_by_bin <- stats::setNames(
    lapply(
      bins,
      function(x) rep(NA_real_, B)
    ),
    bins
  )
  
  if (return_folds) {
    folds_list <- vector(
      "list",
      B
    )
  }
  
  ## ----------------------------------------------------------
  ## Subsampling iterations
  ## ----------------------------------------------------------
  
  for (b in seq_len(B)) {
    
    ## Sample 80% without replacement within every
    ## CHR × AF-bin stratum
    sampled_indices <- unlist(
      Map(
        function(index_vector, n_take) {
          
          sampled_positions <- sample.int(
            n = length(index_vector),
            size = n_take,
            replace = FALSE
          )
          
          index_vector[
            sampled_positions
          ]
        },
        stratum_indices,
        stratum_sample_sizes
      ),
      use.names = FALSE
    )
    
    sub <- d[
      sampled_indices
    ]
    
    ## Check that all chromosomes represented in the complete
    ## dataset remain represented in this iteration
    if (uniqueN(sub$CHR) != n_chr_total) {
      stop(
        paste0(
          "Not all chromosomes were retained in iteration ",
          b,
          "."
        )
      )
    }
    
    ## Correlation across all sampled variants
    rhos_all[b] <- point_rho(sub)
    
    ## Correlation within each AF bin
    tmp <- sub[
      ,
      .(
        rho = point_rho(.SD),
        N = .N,
        n_chr = uniqueN(CHR)
      ),
      by = AF_bin
    ]
    
    ## Store AF-bin correlations
    for (i in seq_len(nrow(tmp))) {
      
      bin_name <- as.character(
        tmp$AF_bin[i]
      )
      
      rhos_by_bin[[bin_name]][b] <- tmp$rho[i]
    }
    
    ## Store fold-level output
    if (return_folds) {
      
      fold_all <- data.table::data.table(
        fold = b,
        AF_bin = "All",
        rho = rhos_all[b],
        N = nrow(sub),
        n_chr = uniqueN(sub$CHR)
      )
      
      fold_bins <- tmp[
        ,
        .(
          fold = b,
          AF_bin,
          rho,
          N,
          n_chr
        )
      ]
      
      folds_list[[b]] <- data.table::rbindlist(
        list(
          fold_all,
          fold_bins
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }
  }
  
  ## ----------------------------------------------------------
  ## Summarize subsampling distributions
  ## ----------------------------------------------------------
  
  summarise_subsampling <- function(
    rhos,
    rho_hat
  ) {
    
    rhos <- rhos[
      is.finite(rhos)
    ]
    
    B_valid <- length(rhos)
    
    if (B_valid == 0L) {
      return(
        tibble::tibble(
          rho = rho_hat,
          p_boot = NA_real_,
          prop_positive = NA_real_,
          prop_negative = NA_real_,
          ci_lo = NA_real_,
          ci_hi = NA_real_
        )
      )
    }
    
    ## Directional subsampling rate.
    ## This is not a conventional independent-sample P value.
    p_boot <- if (alternative == "greater") {
      mean(rhos <= 0)
    } else {
      mean(rhos >= 0)
    }
    
    p_boot <- max(
      p_boot,
      1 / (B_valid + 1)
    )
    
    tibble::tibble(
      rho = rho_hat,
      
      p_boot = p_boot,
      
      prop_positive = mean(
        rhos > 0
      ),
      
      prop_negative = mean(
        rhos < 0
      ),
      
      ci_lo = unname(
        stats::quantile(
          rhos,
          0.025,
          na.rm = TRUE
        )
      ),
      
      ci_hi = unname(
        stats::quantile(
          rhos,
          0.975,
          na.rm = TRUE
        )
      )
    )
  }
  
  ## All variants
  out_all <- summarise_subsampling(
    rhos = rhos_all,
    rho_hat = rho_all_hat
  ) %>%
    dplyr::mutate(
      AF_bin = "All",
      N = nrow(d),
      n_chr = n_chr_total
    )
  
  ## Each AF bin
  out_bins <- lapply(
    names(rhos_by_bin),
    function(bin) {
      
      rho_hat <- rho_bin_hat[
        AF_bin == bin,
        rho
      ]
      
      N_hat <- rho_bin_hat[
        AF_bin == bin,
        N
      ]
      
      n_chr_hat <- rho_bin_hat[
        AF_bin == bin,
        n_chr
      ]
      
      summarise_subsampling(
        rhos = rhos_by_bin[[bin]],
        rho_hat = rho_hat
      ) %>%
        dplyr::mutate(
          AF_bin = bin,
          N = N_hat,
          n_chr = n_chr_hat
        )
    }
  ) %>%
    dplyr::bind_rows()
  
  summary_tbl <- dplyr::bind_rows(
    out_all,
    out_bins
  ) %>%
    dplyr::mutate(
      interval_excludes_zero =
        ci_lo > 0 |
        ci_hi < 0,
      
      direction_label = dplyr::case_when(
        ci_lo > 0 ~ "Robust positive",
        ci_hi < 0 ~ "Robust negative",
        TRUE ~ "Interval includes zero"
      )
    ) %>%
    dplyr::select(
      N,
      n_chr,
      rho,
      p_boot,
      prop_positive,
      prop_negative,
      ci_lo,
      ci_hi,
      interval_excludes_zero,
      direction_label,
      AF_bin
    )
  
  ## ----------------------------------------------------------
  ## Return results
  ## ----------------------------------------------------------
  
  if (return_folds) {
    
    folds_tbl <- data.table::rbindlist(
      folds_list,
      use.names = TRUE,
      fill = TRUE
    )
    
    return(
      list(
        summary = summary_tbl,
        folds = folds_tbl
      )
    )
    
  } else {
    
    return(
      summary_tbl
    )
  }
}


## ============================================================
## 4) Run 2,000 site-subsampling iterations
## ============================================================

res_sitesub <- boot_by_chr_maf_site_spearman(
  df = variants_annotated_clean,
  B = 2000,
  seed = 999,
  alternative = "greater",
  prop_site = 0.8,
  return_folds = TRUE
)

res_sitesub$summary
res_sitesub$folds


data.table::fwrite(
  res_sitesub$summary,
  file =
    "RBA_raw_SiteSubsample_80pctWithinChrxAF_2000fold_OverAllSummary.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

## ============================================================
## 5) Prepare raw fold-level output
## ============================================================

tbl_boot <- data.table::as.data.table(
  res_sitesub$folds
)

tbl_boot[
  ,
  model := "40b_avgRC"
]

data.table::setnames(
  tbl_boot,
  old = c(
    "fold",
    "AF_bin",
    "N",
    "rho"
  ),
  new = c(
    "bootstrap_id",
    "Region",
    "N_variants",
    "rho"
  )
)

tbl_boot <- tbl_boot[
  ,
  .(
    model,
    Region,
    bootstrap_id,
    N_variants,
    rho,
    n_chr
  )
]

data.table::fwrite(
  tbl_boot,
  file =
    "RBA_raw_SiteSubsample_80pctWithinChrxAF_2000fold.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


## ============================================================
## 6) Summarize site-subsampling distributions
## ============================================================

ci_tbl <- tbl_boot[
  ,
  {
    r <- rho[
      is.finite(rho)
    ]
    
    Bn <- length(r)
    
    q2.5 <- if (Bn > 0L) {
      as.numeric(
        stats::quantile(
          r,
          0.025,
          na.rm = TRUE
        )
      )
    } else {
      NA_real_
    }
    
    q97.5 <- if (Bn > 0L) {
      as.numeric(
        stats::quantile(
          r,
          0.975,
          na.rm = TRUE
        )
      )
    } else {
      NA_real_
    }
    
    list(
      B = Bn,
      
      N_each = as.integer(
        round(
          median(
            N_variants,
            na.rm = TRUE
          )
        )
      ),
      
      n_chr_each = as.integer(
        round(
          median(
            n_chr,
            na.rm = TRUE
          )
        )
      ),
      
      rho_mean = if (Bn > 0L) {
        mean(r)
      } else {
        NA_real_
      },
      
      rho_median = if (Bn > 0L) {
        median(r)
      } else {
        NA_real_
      },
      
      rho_sd = if (Bn > 1L) {
        stats::sd(r)
      } else {
        NA_real_
      },
      
      ## Monte Carlo uncertainty of the subsampling mean only
      rho_mc_sem = if (Bn > 1L) {
        stats::sd(r) / sqrt(Bn)
      } else {
        NA_real_
      },
      
      rho_q2.5 = q2.5,
      rho_q97.5 = q97.5,
      
      prop_positive = if (Bn > 0L) {
        mean(r > 0)
      } else {
        NA_real_
      },
      
      prop_negative = if (Bn > 0L) {
        mean(r < 0)
      } else {
        NA_real_
      },
      
      interval_excludes_zero =
        !is.na(q2.5) &&
        !is.na(q97.5) &&
        (
          q2.5 > 0 ||
            q97.5 < 0
        ),
      
      direction_label = if (
        !is.na(q2.5) &&
        q2.5 > 0
      ) {
        "Robust positive"
      } else if (
        !is.na(q97.5) &&
        q97.5 < 0
      ) {
        "Robust negative"
      } else {
        "Interval includes zero"
      }
    )
  },
  by = .(
    model,
    Region
  )
]


## ============================================================
## 7) Write summary output
## ============================================================

data.table::fwrite(
  ci_tbl,
  file =
    "RBA_raw_SiteSubsample_80pctWithinChrxAF_2000fold.summary.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

ci_tbl
