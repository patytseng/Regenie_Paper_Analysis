# Load packages bigsnpr and bigstatsr
library(bigsnpr)

# Load bigSNP objects in R session
print('Loading bigSNP objects...')
for (chr in 1:22) {

  obj.bigSNP <- snp_attach(paste0("/mnt/stsi/stsi1/ptseng/UKBB_Resources/ldpred2/chr",chr,"/temp.rds"))
  #str(obj.bigSNP, max.level = 2, strict.width = "cut")
  
  if (chr == 1) {
    G   <- c(obj.bigSNP$genotypes)
  } else {
    G   <- append(G,obj.bigSNP$genotypes)
  }
}
print('DONE')

# Get map of ld matrix
map <- readRDS("/mnt/stsi/stsi1/ptseng/LD_ref_European_privefl/map.rds")
map$chromosome <- as.numeric(map$chr)

# Get aliases for useful slots
CHR <- map$chromosome
NCORES <- nb_cores()


# Set individuals
set.seed(1)
ind.val <- sample(nrow(G[[1]]), 400000)
ind.test <- setdiff(rows_along(G[[1]]), ind.val)


# Read external summary statistics
print('Loading summary statistics...')
sumstats <- bigreadr::fread2("/gpfs/home/ptseng/Torkamani_Projects/20210712_Analysis-RegeniePaper/exBTs_regenie_phenoCol1_SPA_CADComp.regenie")
#str(sumstats)

names(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "a1_freq", "info", "n_eff", "test", "beta", "beta_se", "chi_sq", "log_p")
#str(obj.bigSNP$map)
map <- map[(1:4)]
names(map) <- c("chr", "pos", "a0", "a1")
#str(map)
#str(sumstats)
sumstats <- snp_match(sumstats, map)
print('DONE')


# Set Temporary File
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# Load Correlation Info from File
print('Loading correlation info...')
for (chr in 1:22) {

  # print(chr)

  ## indices in 'sumstats'
  print('ind-1')
  ind.chr1 <- which(sumstats$chr == chr)
  str(ind.chr1)
  ## indices in 'G'
  print('ind-2')
  ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr1]
  str(ind.chr2)
  ## indices in 'corr'
  print('ind-3')
  ind.chr3 <- match(ind.chr2, which(CHR == chr))
  str(ind.chr3)
    
  print('readRDS')
  corr0 <- readRDS(paste0("/mnt/stsi/stsi1/ptseng/LD_ref_European_privefl/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    df_beta <- sumstats[ind.chr1, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    df_beta <- rbind(df_beta, sumstats[ind.chr1, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
print('DONE')

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
