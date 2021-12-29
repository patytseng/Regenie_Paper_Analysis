jobfile="submit_gctb_1000G_eur.job"
out="gctb_out_r_1000G_eur_sparse"
mkdir ${out}
cd ${out}
mkdir job_reports_gctb_r
echo "

#!/bin/bash
#PBS -N gctb
#PBS -l select=1:ncpus=1:mem=5GB,walltime=00:59:00
#PBS -A UQ-IMB-CNSG
#PBS -o job_reports_gctb_r/
#PBS -e job_reports_gctb_r/
#PBS -J 1-20

out=${out}
" > $jobfile

echo '
cd $PBS_O_WORKDIR
jobno=${PBS_ARRAY_INDEX}
module load intel/2018.3.051 
gctb        --sbayes R \
            --ldm ../../ldm/sparse/chr22/1000G_eur_chr22.ldm.sparse \
            --pi 0.95,0.02,0.02,0.01 \
            --gamma 0.0,0.01,0.1,1 \
            --gwas-summary ../../ma/sim_${jobno}.ma  \
            --out sim_${jobno} \
            --chain-length 10000 --burn-in 2000 --out-freq 10 2>&1 | tee "sim_${jobno}.log"
' >> $jobfile                                                                                        

qsub $jobfile                                                                                      
