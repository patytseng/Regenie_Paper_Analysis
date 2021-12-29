library(bigsnpr)
library(bigreadr)

NCORES <- nb_cores()

## Information for the variants provided in the LD reference
map_ldref <- readRDS("/mnt/stsi/stsi1/ptseng/LD_ref_European_privefl/map.rds")

# Read external summary statistics
print('Loading summary statistics...')
sumstats <- bigreadr::fread2("/gpfs/home/ptseng/Torkamani_Projects/20210712_Analysis-RegeniePaper/exBTs_regenie_phenoCol1_SPA_CADComp.regenie")
names(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "a1_freq", "info", "n_eff", "test", "beta", "beta_se", "chi_sq", "log_p")
#str(sumstats)
print('Done')

info_snp <- snp_match(sumstats, map_ldref)
# 11,792,542 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,054,287 variants have been matched; 0 were flipped and 0 were reversed.
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))

sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))

is_bad <-
  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

library(ggplot2)
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

df_beta <- info_snp[!is_bad, ]

tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {

  cat(chr, ".. ", sep = "")

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

  corr_chr <- readRDS(paste0("/mnt/stsi/stsi1/ptseng/LD_ref_European_privefl/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))
#        int     int_se         h2      h2_se
# 1.09045635 0.01332333 0.12732974 0.01104159
h2_est <- ldsc[["h2"]]

# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, 
                               h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES,
                               sparse = TRUE)  # 13 hours
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
                    
write.csv(beta_auto,'ldpred2_output',row.names=FALSE)
write.csv(df_beta,'ldpred2_statsfile',row.names=FALSE)