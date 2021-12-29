#!/bin/bash

ceildiv(){ echo $((($1+$2-1)/$2)); }

for chr in {22..22}
do

cd full_mcpu/chr$chr

bfile="../../../data/1000G_eur_chr22"
m=`awk -v i=$chr '{if($1==i) print $0}' $bfile.bim | wc -l | cut -d' ' -f1`
k=5000

njob=`ceildiv m k`

echo "
#!/bin/bash

#PBS -N mrg_ch$chr
#PBS -A UQ-IMB-CNSG
#PBS -o stdout
#PBS -e stderr
#PBS -l nodes=1:ppn=1,pmem=40GB,walltime=24:00:00
#PBS -S /bin/bash

chr=$chr

k=$k

filename=1000G_eur_chr${chr}
" > merge.job

echo '
cd $PBS_O_WORKDIR

out=$filename

> $out.mldmlist
' >> merge.job

echo "
for i in {1..$njob}
do
" >> merge.job

echo '
echo "${PWD}/$out.snp$((k*(i-1)+1))-$((k*i)).ldm.full" >> $out.mldmlist

done
module load intel/2018.3.051
gctb --mldm $out.mldmlist --make-full-ldm --out $out > $out.log 2>&1
' >> merge.job

#qsub merge.job

cd ../../
done
