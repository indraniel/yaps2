#!/bin/bash

# Inputs:
#    datadir ($1) : the top level directory containing the per chromosome variant eval stats
#                   (e.g. $basedir/7-gatk-variant-eval)
#    outdir  ($2) : the output directory that will contain the merged variant stats and plots
#                   (e.g. $basedir/7.1-gatk-variant-eval-summary)


# setup env
BIO_1662=/gscmnt/gc2802/halllab/idas/jira/BIO-1662
cd $BIO_1662
source .env

echo "# ------- merge GATK variant eval results -------";
datadir=$1
outdir=$2
merged_stats=${outdir}/merged.dat

python ${BIO_1662}/bin/merge-gatk-vcfeval-outputs.py \
    --out ${merged_stats} \
    --skip 'MetricsCollection' \
    $(find ${datadir} -name "*eval.out")

echo "# ------- plot merged GATK variant eval results -------";
pdf=${outdir}/gatk-variant-eval-stats.pdf

Rscript ${BIO_1662}/bin/diagnostic-plots-annotation-review.r --input=${merged_stats} --output=${pdf}
echo "created ${pdf}"

echo "# ------- plot merged GATK variant eval results (minus controls) -------";
pdf=${outdir}/gatk-variant-eval-stats-minus-control-samples.pdf

Rscript ${BIO_1662}/bin/diagnostic-plots-annotation-review.r \
    --input=${merged_stats} \
    --prune="H_HF-CHM1htert-US-5A,H_IJ-NA12878-NA12878_K5,H_PY-CHM13-CHM13h" \
    --output=${pdf}

echo "created ${pdf}"
