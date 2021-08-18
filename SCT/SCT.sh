for chr in {1..22};
    do sbatch \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=20 \
        --time=7-00:00:00 \
        --mem=200GB \
        --partition=highmem \
        --export=chr=${chr} \
        --job-name=SCT_chr${chr} \
        --out=SCT_chr${chr}_jobID%j \
        SCT.qsub;
done