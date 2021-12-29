for chr in {1..22}; \
    do \
    sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=16 \
    --time=7-00:00:00 \
    --mem=60GB \
    --export=chr=${chr} \
    --job-name=plink_pt_chr${chr} \
    --out=SCT_chr${chr}_jobID%j \
    SCT.qsub; \
    done