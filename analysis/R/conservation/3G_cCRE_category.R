library(readr)

data=read_tsv("Evo2_6Models_Inference_result.6475639.UCSC_Fullanno.LCR_addMAF_Tier.June10.txt.gz")

cols <- c(
  "rsid",
  "BP",
  "Ref",
  "Alt",
  "Alt_Freq",
  "40b_avgRC_delta_score",
  "Func.ensGene",
  "Gene.ensGene",
  "GeneDetail.ensGene",
  "category",
  "AAChange.ensGene",
  "phastCons100way",
  "phyloP100way",
  "ExonicFunc.ensGene",
  'repName'
)


data_clean=data[,cols]
data_clean=data_clean[is.na(data_clean$repName),]

table(data_clean$category)

data_clean=data_clean[data_clean$Ref %in% c("A","T","C","G") & data_clean$Alt %in% c("A","T","C","G"),]

### All by coding noncoding

# Example:
# df$score  -> numeric column
# df$group  -> factor column

# Make sure grouping column is factor
data_clean$group <- as.factor(data_clean$category)

### Test for value 

library(dplyr)

dat <- data_clean %>%
  mutate(
    score = abs(`40b_avgRC_delta_score`),
    cat   = as.factor(category)
  ) %>%
  filter(!is.na(score), !is.na(cat))

ref <- "Intergenic"
cats_to_test <- setdiff(levels(dat$cat), ref)

res <- lapply(cats_to_test, function(g) {
  
  d2 <- dat %>% filter(cat %in% c(ref, g))
  
  # Wilcoxon
  wt <- wilcox.test(score ~ cat, data = d2, exact = FALSE)
  
  # Summary statistics
  mean_ref   <- mean(d2$score[d2$cat == ref])
  mean_grp   <- mean(d2$score[d2$cat == g])
  
  median_ref <- median(d2$score[d2$cat == ref])
  median_grp <- median(d2$score[d2$cat == g])
  
  data.frame(
    group = g,
    n_ref = sum(d2$cat == ref),
    n_grp = sum(d2$cat == g),
    
    mean_ref   = mean_ref,
    mean_group = mean_grp,
    mean_diff  = mean_grp - mean_ref,
    
    median_ref   = median_ref,
    median_group = median_grp,
    median_diff  = median_grp - median_ref,
    
    p_raw = wt$p.value
  )
  
}) |> bind_rows()

# Adjust p-values (fdr)
res$p_adj <- p.adjust(res$p_raw, method = "fdr")

# Order by adjusted p
res <- res[order(res$p_adj), ]

res


summary_all <- dat %>%
  group_by(cat) %>%
  summarise(
    n = n(),
    mean_score   = mean(score),
    median_score = median(score),
    sd_score     = sd(score),
    .groups = "drop"
  )

summary_all
