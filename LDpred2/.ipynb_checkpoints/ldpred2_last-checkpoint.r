# Load packages bigsnpr and bigstatsr
library(bigsnpr)

# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("/mnt/stsi/stsi1/ptseng/UKBB_Resources/ukb41999_imp_v3.rds")
obj.bigSNP$map$chromosome <- as.numeric(obj.bigSNP$map$chromosome)
# See how the file looks like
#str(obj.bigSNP, max.level = 2, strict.width = "cut")

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
NCORES <- nb_cores()


# Set individuals
set.seed(1)
ind.val <- sample(nrow(G), 400000)
ind.test <- setdiff(rows_along(G), ind.val)


# Read external summary statistics
sumstats <- bigreadr::fread2("/gpfs/home/ptseng/Torkamani_Projects/20210602_regenie-UKBB/regenie-ukbb_impute4_step2_f.50.0.0.regenie")
#str(sumstats)

names(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "a1_freq", "info", "n_eff", "test", "beta", "beta_se", "chi_sq", "log_p")
#str(obj.bigSNP$map)
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
#str(map)
sumstats <- snp_match(sumstats, map)


# Set Temporary File
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# Load Correlation Info from File
for (chr in 1:22) {

  # print(chr)

  ## indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  ## indices in 'G'
  ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
  ## indices in 'corr'
  ind.chr3 <- match(ind.chr2, which(CHR == chr))

  corr0 <- readRDS(paste0("/gpfs/home/ptseng/Torkamani_Projects/20210602_regenie-UKBB/corr/chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    df_beta <- sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    df_beta <- rbind(df_beta, sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]

# ldpred2 auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.val,
                         ind.col = df_beta[["_NUM_ID_"]])
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])

# compute predictions for test set
pred_test <- big_prodMat(G, final_beta_auto, ind.row = ind.test,
                         ind.col = df_beta[["_NUM_ID_"]])

# save results
res <- list(pred = setNames(as.data.frame(pred_test), colnames(betas)),
            params = params, auto = multi_auto[keep])
saveRDS(res, '/gpfs/home/ptseng/Torkamani_Projects/20210602_regenie-UKBB/ldpred2_result')
