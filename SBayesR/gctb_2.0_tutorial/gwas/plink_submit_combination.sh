#!/bin/bash
#PBS -N plink_gwas_1kg
#PBS -o stdout
#PBS -e stderr
#PBS -l select=1:ncpus=1:mem=1GB,walltime=00:10:00
#PBS -A UQ-IMB-CNSG
#PBS -J 1-20


cd $PBS_O_WORKDIR
jobno=${PBS_ARRAY_INDEX}

plink --bfile ../data/1000G_eur_chr22 \
      --linear hide-covar \
      --allow-no-sex \
      --pheno ../pheno/phenos_ga1.txt  \
      --mpheno ${jobno} \
      --out sim_${jobno}  2>&1 | tee "sim_${jobno}.log"

