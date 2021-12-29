for chr in {1..22};
    do sbatch \
        --nodes=1 \
        --exclude=nodea[1301-1323,1325-1331],emb[0702-0736],imsb[0501-0532,0601-0632],nodeb[1219-1228] \
        --ntasks=1 \
        --cpus-per-task=16 \
        --time=7-00:00:00 \
        --mem=120GB \
        --partition=shared \
        --export=chr=${chr} \
        --job-name=SBayesR_chr${chr} \
        --out=%x_jobID%j \
        sbayesr.qsub;
done

# emb[0702-0736]
# imsb[0501-0532,0601-0632]
# nodea[1301-1323,1325-1331]