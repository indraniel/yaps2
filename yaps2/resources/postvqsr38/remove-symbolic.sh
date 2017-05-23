#!/bin/bash

set -ueo pipefail

BCFTOOLS=/gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4
TABIX=/gscmnt/gc2802/halllab/idas/software/local/bin/tabix

INVCF=$1
OUTVCF=$2

if [ -e $OUTVCF ]
then
    exit 0
fi

TMPVCF=$OUTVCF.temp
${BCFTOOLS} view -e '%TYPE="other" || ALT="*"' $INVCF --output-type z --output-file $TMPVCF \
    && ${TABIX} -p vcf -f $TMPVCF && mv $TMPVCF.tbi $OUTVCF.tbi && mv $TMPVCF $OUTVCF
