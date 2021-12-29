# Load packages
library(bigsnpr)
library(bigassertr)
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

# Attach the "bigSNP" object in R session
print('Grabbing bigSNP Object...')
obj.bigSNP <- snp_attach(paste0("/mnt/stsi/stsi1/ptseng/UKBB_Resources/ldpred2/chr",chr,"/temp.rds"))
obj.bigSNP$map$chromosome <- as.numeric(obj.bigSNP$map$chromosome)
# See how the file looks like
#str(obj.bigSNP, max.level = 2, strict.width = "cut")
print('Grabbing bigSNP Object... done')

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

POS2 <- snp_asGeneticPos(CHR, POS, dir = "/gpfs/home/ptseng/Torkamani_Projects/20210712_Analysis-RegeniePaper", ncores = 16)

# Generate Correlation Info into File
ind.chr <- which(CHR == chr)

print('Starting to make LD Matrix')

runonce::save_run(
    snp_cor(
        G, ind.row = ind.val, ind.col = ind.chr,
        alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000,
        ncores = NCORES
    ),
    file = paste0("/gpfs/home/ptseng/Torkamani_Projects/20210712_Analysis-RegeniePaper/corr/chr", chr, ".rds")
)