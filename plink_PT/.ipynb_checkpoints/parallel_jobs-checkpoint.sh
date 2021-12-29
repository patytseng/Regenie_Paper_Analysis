for chrom in {1..22}; \
    do \
    echo sbatch \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=16 \
    --mem=120gb \
    --export=path=/mnt/stsi/stsi3/External/UKBB/genomic_data/geneotypes/chr${chrom}/ukb_imp_chr${chrom}_v3.bgen \
    --job-name=plink_pt_chr${chrom} \
    --out=plink_pt_chr${chrom}_%j \
    plink_pt.qsub; \
    done > files_to_merge.txt
