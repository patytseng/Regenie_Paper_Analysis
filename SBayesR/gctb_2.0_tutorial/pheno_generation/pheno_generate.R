# --------------------------------------------
# Generate 20 simulation replicates with 
# using 1000G data with a h2_snp=0.1 and
# 1500 causal variants 2 with expalin 2% and
# 3% of the genetic variance respectively.
# Author: Luke Lloyd-Jones
# Date: 29/01/2019
# --------------------------------------------
# First read the data and save as an RData 
# object for easy data processing
source("read_plink_func.R")
# Read the bed file
onekg_22 <- read.plink("../data/1000G_eur_chr22")
# Read the bim and fam file and add the row and column names
# to the matrix
bim <- read.table("../data/1000G_eur_chr22.bim")
fam <- read.table("../data/1000G_eur_chr22.fam") 
rownames(onekg_22) <- as.character(fam[, 1])
colnames(onekg_22) <- as.character(bim[, 2])
# Cycle over the columns and center and scale under HWE
ninds <- dim(onekg_22)[1]
nmarkers <- dim(onekg_22)[2]
afs <- array(0,dim(onekg_22)[2])
X <- matrix(0, ncol = nmarkers, nrow = ninds)
rownames(X) <- as.character(fam[, 1])
colnames(X) <- as.character(bim[, 2])	
for (j in seq(1, dim(onekg_22)[2]))
{
 print(paste0("Doing marker ", j))
 # Calculate the reference allele frequqnecy
 nas <- which(is.na(onekg_22[, j]))
 if (length(as.numeric(nas)) > 0)
 {
   q <- sum(onekg_22[-nas, j], na.rm = T) / (2 * (ninds- length(as.numeric(nas))))
 } else {
   q <- sum(onekg_22[, j], na.rm = T) / (2 * (ninds))
 }
 afs[j] <- q		
 # Centre and then scale	    
 X[, j] <- (onekg_22[, j] - 2 * q) / sqrt(2 * q * (1 - q))
}
# --------------------------------------------
# Check that the matrix is centered and scaled
# --------------------------------------------
# They should all have an af
summary(afs)
afs
colMeans(X, na.rm = T)
apply(X, 2, function(x) { sd(x, na.rm = T)})
# Write the centered matrix out
save(X, file = "1000G_genos_scaled.Rdata")
# ===============================================
# Phenotype 1
# ===============================================
load("1000G_genos_scaled.Rdata")
nmarkers <- dim(X)[2]
ninds    <- dim(X)[1]
fam <- read.table("../data/1000G_eur_chr22.fam")
phenos <- fam[,1:2]
for (i in seq(1, 20))
{
  print(paste0("Phenotype ", i)) 
  # Sample the 1000 tail variants
  variants   <- sample(seq(1, nmarkers), 1500)         
  variants.o <- variants
  # Use the first three as the big three
  beta1 <- sqrt(0.03)
  beta2 <- sqrt(0.02)
  # Poly tail
  beta.poly.tail <- rnorm(1498, 0, sqrt(0.05/1498))
  # Sum to get G
  g <- X[, c(variants.o)] %*% c(c(beta1, beta2), beta.poly.tail)
  h2=0.1
  y <- g + rnorm(ninds, 0, sqrt(var(g)*(1/h2 - 1)))
  y2 <- (y - mean(y)) / sd(y)
  print(summary(lm(y~X[, c(variants.o[1:10])])))
  # Write out the variances
  vars <- data.frame(var(g), var(y), var(g)/var(y))
  vars.out <- paste0("../pheno/variances_ga1_sim_", i, ".txt")
  colnames(vars) <- c("VAR_G", "VAR_Y", "H2")
  write.table(vars, vars.out, col.names = T, row.names = F, sep = "\t", quote = F)
  # Bind up the phenotypes
  phenos <- cbind(phenos, y2) 
  # Write out the SNP effects
  snps <- colnames(X)[c(variants.o)]
  snps.effs <- cbind(snps, c(c(beta1, beta2), beta.poly.tail))
  snps.effs.out <- paste0("../pheno/snp_effects_ga1_sim_", i, ".txt")
  colnames(snps.effs) <- c("SNP", "EFFECT") 
  write.table(snps.effs, snps.effs.out, col.names = T, row.names = F, sep = "\t", quote = F) 
}
phenos.out       <- paste0("../pheno/phenos_ga1.txt")
write.table(phenos, phenos.out, col.names = F, row.names = F, sep = " ", quote = F)











