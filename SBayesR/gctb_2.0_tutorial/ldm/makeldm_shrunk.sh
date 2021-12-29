#!/bin/bash


ceildiv(){ echo $((($1+$2-1)/$2)); }

for chr in {22..22}
do

mkdir -p shrunk/chr$chr
cd shrunk/chr$chr
mkdir stderr
mkdir stdout
bfile="../../../data/1000G_eur_chr22" 

jobfile="makeldm.job"

echo "
#!/bin/bash

#PBS -N ldm_ch${chr}
#PBS -A UQ-IMB-CNSG
#PBS -o stdout
#PBS -e stderr
#PBS -l select=1:ncpus=1,pmem=10GB,walltime=00:30:00
#PBS -S /bin/bash

chr=$chr

bfile=${bfile}
out=1000G_eur_chr${chr}
" > $jobfile

echo '
cd $PBS_O_WORKDIR

module load intel/2018.3.051  
gctb --bfile $bfile \
     --make-shrunk-ldm \
     --gen-map ../../genetic_map/chr22.OMNI.interpolated_genetic_map \
     --write-ldm-txt \
     --out $out > $out.ldm.shrunk.log 2>&1   
' >> $jobfile

qsub $jobfile

cd ../../
done
