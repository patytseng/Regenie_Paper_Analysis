# Load packages bigsnpr and bigstatsr
library(bigsnpr)
library(stringr)
library(optparse)

# Obtaining chromosome number
option_list = list(
  make_option("--chr", type="integer", default=NULL, help="chromosome numbers")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$chr)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (chromosome number).n", call.=FALSE)
}

chr <- opt$chr

# Grabbing SNP IDs from stats file
print('Grabbing SNP IDs from stats file')
data <- read.csv('exBTs_regenie_phenoCol1_SPA_CADComp.regenie',sep=' ')
data <- data[,c("CHROM","GENPOS","ALLELE0","ALLELE1")]
data["INDEX"] = str_c(data$CHROM,data$GENPOS,data$ALLELE0,data$ALLELE1,sep="_")

# Formatting list
print('Finishing list of SNP IDs')
list_snp_id <- list()
for (val in chr) {
    subset <- data[(data$CHROM == val),]
    list_snp_id <- append(list_snp_id, list(subset$INDEX))
}

# Read from BGEN, it generates .bk and .rds files
print('Reading from BGEN and generating bigSNP object with selected SNPs')
rds <- snp_readBGEN(
  bgenfiles = glue::glue("/mnt/stsi/stsi1/ptseng/UKBB_Resources/ukbb_impute4/bgen_links/ukb_imp_chr{chr}_v3.bgen", chr = chr),
  backingfile = "/scratch_ssd/temp",
  list_snp_id = list_snp_id,
  ncores = 1
)
