#!/bin/bash

set -ueo pipefail

BIO_1662=/gscmnt/gc2802/halllab/idas/jira/BIO-1662
BCFTOOLS=/gscmnt/gc2719/halllab/bin/bcftools
TABIX=/gscmnt/gc2719/halllab/bin/tabix

INVCF=$1
OUTVCF=$2

if [ -a $OUTVCF ]
then
    exit 0
fi

REF=/gscmnt/gc2719/halllab/genomes/human/GRCh37/1kg_phase1/human_g1k_v37.fasta
TMPVCF=$OUTVCF.temp
${BCFTOOLS} view -e '%TYPE="other"' $INVCF --output-type z --output-file $TMPVCF \
    && ${TABIX} -p vcf -f $TMPVCF && mv $TMPVCF.tbi $OUTVCF.tbi && mv $TMPVCF $OUTVCF
