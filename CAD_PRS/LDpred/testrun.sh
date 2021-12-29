for weights in /gpfs/home/ptseng/Torkamani_Projects/20210712_Analysis-RegeniePaper/CAD_PRS/LDpred/reweighted_LDpred_V*.csv;
    do sbatch \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=16 \
        --time=7-00:00:00 \
        --mem=60GB \
        --export=weights=${weights} \
        --job-name=testrun_PRS_validator \
        --out=testrun_PRS_validator_jobID%j \
        testrun.qsub; \
done