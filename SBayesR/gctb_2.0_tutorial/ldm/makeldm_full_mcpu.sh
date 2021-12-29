#!/bin/bash


ceildiv(){ echo $((($1+$2-1)/$2)); }

for chr in {22..22}
do

mkdir -p full_mcpu/chr$chr
cd full_mcpu/chr$chr
mkdir stderr
mkdir stdout
bfile="../../../data/1000G_eur_chr22" 
m=`awk -v i=$chr '{if($1==i) print $0}' $bfile.bim | wc -l | cut -d' ' -f1`
k=5000

njob=`ceildiv m k`

jobfile="makeldm.job"

echo "
#!/bin/bash

#PBS -N ldm_ch${chr}
#PBS -A UQ-IMB-CNSG
#PBS -o stdout
#PBS -e stderr
#PBS -l select=1:ncpus=1,pmem=10GB,walltime=00:10:00
#PBS -S /bin/bash
#PBS -J 1-$njob

chr=$chr

k=$k
bfile=${bfile}
out=1000G_eur_chr${chr}
" > $jobfile

echo '
cd $PBS_O_WORKDIR

i=${PBS_ARRAY_INDEX}
module load intel/2018.3.051  
gctb --bfile $bfile --make-full-ldm --snp $((k*(i-1)+1))-$((k*i)) --out $out > $out.snp$((k*(i-1)+1))-$((k*i)).ldm.full.log 2>&1   
' >> $jobfile

#qsub $jobfile

cd ../../
done
