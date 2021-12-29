# --------------------------------------------
# Generate 10 simulation replicates with 
# three variants explaing 10, 5, and 2% of the 
# PV and a ply tail with 1000 variants
# Date: 02/10/2018
# --------------------------------------------
# First read the data and save as an RData 
# object for easy data processing
source("read_plink_func.R")
# Read the bed file
ukb_2122 <- read.plink("../data/ukbEURu_hm3_v3_10k_chr2122_valid")
# Read the bim and fam file and add the row and column names
# to the matrix
bim <- read.table("../data/ukbEURu_hm3_v3_10k_chr2122_valid.bim")
fam <- read.table("../data/ukbEURu_hm3_v3_10k_chr2122_valid.fam") 
rownames(ukb_2122) <- as.character(fam[, 1])
colnames(ukb_2122) <- as.character(bim[, 2])
# Cycle over the columns and center and scale under HWE
ninds <- dim(ukb_2122)[1]
nmarkers <- dim(ukb_2122)[2]
afs <- array(0,dim(ukb_2122)[2])
X <- matrix(0, ncol = nmarkers, nrow = ninds)	
for (j in seq(1, dim(ukb_2122)[2]))
{
 print(paste0("Doing marker ", j))
 # Calculate the reference allele frequqnecy
 nas <- which(is.na(ukb_2122[, j]))
 if (length(as.numeric(nas)) > 0)
 {
   q <- sum(ukb_2122[-nas, j], na.rm = T) / (2 * (ninds- length(as.numeric(nas))))
 } else {
   q <- sum(ukb_2122[, j], na.rm = T) / (2 * (ninds))
 }
 afs[j] <- q		
 # Centre and then scale	    
 X[, j] <- (ukb_2122[, j] - 2 * q) / sqrt(2 * q * (1 - q))
}
# --------------------------------------------
# Check that the matrix is centered and scaled
# --------------------------------------------
# They should all have an af
summary(afs)
rownames(X) <- as.character(fam[, 1])
colnames(X) <- as.character(bim[, 2])
colMeans(X, na.rm = T)
apply(X, 2, function(x) { sd(x, na.rm = T)})
# ------------------------------------------------
# Impute the mean for prediction generation
# ------------------------------------------------
ninds <- dim(X)[1]
nmarkers <- dim(X)[2]
X[is.na(X)] <- 0 
save(X, file = "ukb_eu_hm3_chr21_22_scaled_imp_validation.Rdata")
# ===============================================
# Phenotype 1
# ===============================================
load("../pheno1/ukb_eu_hm3_chr21_22_scaled_imp_validation.Rdata")
fam <- read.table("../data/ukbEURu_hm3_v3_10k_chr2122_valid.fam")
ninds <- dim(X)[1]
phenos <- fam[,1:5]
for (i in seq(1, 10))
{
  print(paste0("Phenotype ", i)) 
  # Read in the effects
  effs <- read.table(paste0("ga1/snp_effects_ga1_sim_", i, ".txt"), header = T)
  # Make the breeding values then the phenotype
  g <- X[, as.character(effs$SNP)] %*% effs$EFFECT
  h2=0.1 
  y <- g + rnorm(ninds, 0, sqrt(var(g)*(1/h2 - 1)))
  y2 <- (y - mean(y)) / sd(y) 
  # Bind up the phenotypes
  phenos <- cbind(phenos, y2) 
}
phenos.bolt <- phenos[, c(1, 2, seq(6,dim(phenos)[2]))]
phenos.bolt.out  <- paste0("phenos_ga1_bolt_validation.txt")
phenos.out       <- paste0("phenos_ga1_validation.txt")
colnames(phenos.bolt) <- c("FID", "IID", paste0("SIM_", seq(1, 10)))
write.table(phenos.bolt, phenos.bolt.out, col.names = T, row.names = F, sep = " ", quote = F)
write.table(phenos, phenos.out, col.names = F, row.names = F, sep = " ", quote = F)
# ===============================================
# Phenotype 2 - R model 
# ===============================================
phenos <- fam[,1:5]
for (i in seq(1, 10))
{ 
  print(paste0("Phenotype ", i))                                                                                    
  # Read in the effects
  effs <- read.table(paste0("ga2/snp_effects_ga2_sim_", i, ".txt"), header = T)                                     
  # Make the breeding values then the phenotype
  g <- X[, as.character(effs$SNP)] %*% effs$EFFECT
  h2=0.1 
  y <- g + rnorm(ninds, 0, sqrt(var(g)*(1/h2 - 1)))
  y2 <- (y - mean(y)) / sd(y)                                                                                       
  # Bind up the phenotypes                                                                                          
  phenos <- cbind(phenos, y2)                                                                                       
}                                                                                                                   
phenos.bolt <- phenos[, c(1, 2, seq(6,dim(phenos)[2]))]                                                             
phenos.bolt.out  <- paste0("phenos_ga2_bolt_validation.txt")                                                            phenos.out       <- paste0("phenos_ga2_validation.txt")                                                                 colnames(phenos.bolt) <- c("FID", "IID", paste0("SIM_", seq(1, 10)))
write.table(phenos.bolt, phenos.bolt.out, col.names = T, row.names = F, sep = " ", quote = F)
write.table(phenos, phenos.out, col.names = F, row.names = F, sep = " ", quote = F)
# ===============================================
# Phenotype 3 - BLUP model 
# ===============================================
phenos <- fam[,1:5]
for (i in seq(1, 10))
{
  print(paste0("Phenotype ", i))                                                                                    
  # Read in the effects                                                                                             
  effs <- read.table(paste0("ga3/snp_effects_ga3_sim_", i, ".txt"), header = T) 
  # Make the breeding values then the phenotype
  g <- X[, as.character(effs$SNP)] %*% effs$EFFECT
  h2=0.1 
  y <- g + rnorm(ninds, 0, sqrt(var(g)*(1/h2 - 1)))
  y2 <- (y - mean(y)) / sd(y)
  # Bind up the phenotypes                                                                                          
  phenos <- cbind(phenos, y2)
}
phenos.bolt <- phenos[, c(1, 2, seq(6,dim(phenos)[2]))]
phenos.bolt.out  <- paste0("phenos_ga3_bolt_validation.txt")                                                            phenos.out       <- paste0("phenos_ga3_validation.txt")
colnames(phenos.bolt) <- c("FID", "IID", paste0("SIM_", seq(1, 10)))
write.table(phenos.bolt, phenos.bolt.out, col.names = T, row.names = F, sep = " ", quote = F)
write.table(phenos, phenos.out, col.names = F, row.names = F, sep = " ", quote = F)











