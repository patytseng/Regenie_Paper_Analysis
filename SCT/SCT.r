# Load packages
library(bigsnpr)
library(optparse)

# Obtaining chromosome number
option_list = list(
    make_option("--chr", type="integer", default=NULL, help="chromosome numbers")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$chr)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (chromosome number)", call.=FALSE)
}

chr <- opt$chr

cat(paste0('\nCurrently running SCT on chromosome ',chr,'\n\n'))

# Reading bigSNP
obj.bigSNP <- snp_attach("/scratch_ssd/temp.rds")
obj.bigSNP$map$chromosome <- as.numeric(obj.bigSNP$map$chromosome)

# Get aliases for useful slots
G <- obj.bigSNP$genotypes
map <- obj.bigSNP$map
CHR <- map$chromosome
POS <- map$physical.pos
NCORES <- nb_cores()

# Read external phenotype file
phenotype <- read.csv('/mnt/stsi/stsi1/ptseng/UKBB_Resources/phenos/CAD.csv')

# Matching individuals between genotypes and phenotypes
genotype_ID <- read.csv('/mnt/stsi/stsi1/ptseng/UKBB_Resources/ID_list/impute4.txt')
names(genotype_ID) <- 'ID'
phenotype$use <- TRUE
matched_phenotype <- merge(genotype_ID,phenotype,by = 'ID',all.x = TRUE,sort=FALSE)
matched_phenotype <- merge(genotype_ID,matched_phenotype,by = 'ID',sort=FALSE)

# Get aliases for useful slots
y <- as.vector(na.omit(matched_phenotype$CAD_Composite))
IND <- which(matched_phenotype$use)

# Read external summary statistics file
sumstats <- bigreadr::fread2("/gpfs/home/ptseng/Torkamani_Projects/20210712_Analysis-RegeniePaper/exBTs_regenie_phenoCol1_SPA_CADComp.regenie")

# Matching summary statistics with bigSNP map
names(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "a1_freq", "info", "n_eff", "test", "beta", "beta_se", "chi_sq", "log_p")
minimap <- map[,-(2:3)]
names(minimap) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats[sumstats$chr == chr,], minimap)

# beta and lpval need to have the same length as ncol(G), CHR and POS
beta <- rep(NA, ncol(G))
beta[info_snp$'_NUM_ID_'] <- info_snp$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp$'_NUM_ID_'] <- info_snp$log_p

cat('\nBegninning Clumping\n\n')
# Clumping
all_keep <- snp_grid_clumping(
    G = G, 
    infos.chr = CHR, 
    infos.pos = POS, 
    lpS = lpval, 
    ind.row = IND,
    exclude = which(is.na(lpval)),
    ncores = NCORES
)
attr(all_keep, "grid")

cat('\nBegninning Thresholding')
# Thresholding
multi_PRS <- snp_grid_PRS(
    G = G, 
    all_keep = all_keep, 
    betas = beta, 
    lpS = lpval, 
    n_thr_lpS = 50, 
    ind.row = IND,
    backingfile = "/scratch_ssd/thresholding_temp", 
    ncores = NCORES
)

cat('\nBegninning Stacking')
# Stacking
final_mod <- snp_grid_stacking(
    multi_PRS = multi_PRS, 
    y.train = y, 
    K = 10, 
    ncores = NCORES
)

cat('\nProcessing Output')
# Making weights table
output <- data.frame(info_snp$chr, info_snp$pos, info_snp$a0, info_snp$a1, na.omit(final_mod$beta.G))
names(output) <- c('chr','pos','a0','a1','beta')

cat('\nWriting Output')
# Saving weights table
write.csv(output,paste0('chr',chr,'_weights.csv'),row.names = FALSE)
