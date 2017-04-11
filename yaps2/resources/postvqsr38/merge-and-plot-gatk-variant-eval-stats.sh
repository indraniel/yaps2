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

if [ ! -d "$outdir" ]; then
    echo "Creating output directory : $outdir"
    mkdir -p $outdir
fi

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
    --prune="H_HF-CHM1htert-US-5A,H_IJ-HG00512-HG00512_1,H_IJ-HG00513-HG00513_1,H_IJ-HG00514-HG00514_1,H_IJ-HG00731-HG00731_2,H_IJ-HG00732-HG00732_1,H_IJ-HG00733-HG00733_2,H_IJ-HG01350-HG01350_1,H_IJ-HG01351-HG01351_1,H_IJ-HG01352-HG01352_1,H_IJ-HG02059-HG02059_1,H_IJ-HG02060-HG02060_1,H_IJ-HG02061-HG02061_1,H_IJ-HG02816-HG02816_1,H_IJ-HG02817-HG02817_1,H_IJ-HG02818-HG02818_1,H_IJ-NA12878-NA12878_K10,H_IJ-NA12891-NA12891_D2,H_IJ-NA12892-NA12892_E1,H_IJ-NA19238-NA19238_D3,H_IJ-NA19239-NA19239_B9,H_IJ-NA19240-NA19240_F1,H_IJ-NA24143-NA24143_B2,H_IJ-NA24149-NA24149_B1,H_IJ-NA24385-NA24385_3,H_PY-CHM13-CHM13h" \
    --output=${pdf}

echo "created ${pdf}"
