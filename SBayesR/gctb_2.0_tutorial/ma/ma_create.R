library(dplyr)
args = commandArgs(trailingOnly=TRUE)
for (i in seq(1, args[3]))
{
  print(i)
  x.in <- paste0(args[1], ".bim")
  y.in <- paste0(args[1], ".freq")
  z.in <- paste0(args[2], "sim_", i, ".assoc.linear")
  x <- read.table(x.in)
  y <- read.table(y.in)
  z <- read.table(z.in, header = T)
  x.y   <- inner_join(x, y, by = c("V2" = "V1"))
  x.y.z <- inner_join(x.y, z,  by = c("V2" = "SNP"))
  x.y.z$SE <- x.y.z$BETA/x.y.z$STAT
  x.y.z.2 <- x.y.z[, c("V2", "V5", "V6", "V3.y", "BETA", "SE", "P", "NMISS")]
  colnames(x.y.z.2) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
  write.table(x.y.z.2, paste0(args[4], "sim_", i, ".ma"), 
              col.names = T, row.names = F, sep = " ", quote = F)
}

