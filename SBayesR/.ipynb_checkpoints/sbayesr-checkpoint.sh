for chr in {1..22};
    do sbatch \
        --nodes=1 \
        --nodelist='nodeb[0201-1228]' \
        --ntasks=1 \
        --cpus-per-task=16 \
        --time=7-00:00:00 \
        --mem=100GB \
        --partition=shared \
        --export=chr=${chr} \
        --job-name=SBayesR_chr${chr} \
        --out=%x_jobID%j \
        sbayesr.qsub;
done