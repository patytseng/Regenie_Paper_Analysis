# -------------------------------------------
# An R code script to read plink binary files
# -------------------------------------------

read.plink <- function(root)
{
  bed.file = paste(root, ".bed", sep = "")
  bed.file.size = file.info(bed.file)$size
  sample.size = dim(read.table(paste(root, ".fam", sep = "")))[1]
  snp.size = ceiling(sample.size/4)
  n.snps = round((bed.file.size - 3)/snp.size)
  bin.connection = file(bed.file, "rb")
  test.bytes = readBin(bin.connection, what = "raw", n = 3)
  if (!identical(as.character(test.bytes), c("6c", "1b", "01"))) {
    stop("BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file")
  }
  genotypes = matrix(ncol = n.snps, nrow = sample.size)
  for (i in 1:n.snps) {
    r.bin.snp = readBin(bin.connection, what = "raw", n = snp.size)
    bin.snp = matrix(as.numeric(rawToBits(r.bin.snp)), ncol = 2,
                     byrow = T)[1:sample.size, ]
    genotypes[, i] = bin.snp[, 1] + bin.snp[, 2] - 10 * ((bin.snp[, 1] == 1) & (bin.snp[, 2] == 0))
  }
  genotypes[genotypes == -9] = NA
  #snp.names = read.table(paste(root, ".bim", sep = ""))[, 2]
  #colnames(genotypes) = snp.names
  close(bin.connection)
  return(genotypes)
}
