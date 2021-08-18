for chr in {1..22};
    do sbatch \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=16 \
        --time=7-00:00:00 \
        --mem=100GB \
        --partition=shared \
        --export=chr=${chr} \
        --job-name=PRScs_chr${chr} \
        --out=PRScs_chr${chr}_jobID%j \
        PRScs.qsub;
done