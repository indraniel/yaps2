#!/bin/bash

set -ueo pipefail

BCFTOOLS=/gscmnt/gc2719/halllab/bin/bcftools
TABIX=/gscmnt/gc2719/halllab/bin/tabix

INVCF=$1
OUTVCF=$2

if [ -e $OUTVCF ]
then
    exit 0
fi

TMPVCF=$OUTVCF.temp
zcat $INVCF | awk '$5!="<*:DEL>" && $5!="*"' | bgzip -c > $TMPVCF \
    && ${TABIX} -p vcf -f $TMPVCF && mv $TMPVCF.tbi $OUTVCF.tbi && mv $TMPVCF $OUTVCF
