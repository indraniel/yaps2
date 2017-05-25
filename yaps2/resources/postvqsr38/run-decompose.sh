#!/bin/bash

# modified from @aregier
#  ~aregier/scratch/dlarson-gatk-scripts/run_decompose.sh

BIO_1662=/gscmnt/gc2802/halllab/idas/jira/BIO-1662
VT=${BIO_1662}/vendor/local/bin/vt-0.5
TABIX=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin/tabix

INVCF=$1
OUTVCF=$2
CHROM=$3

if [ -a $OUTVCF ]
then
    exit 0
fi

REF=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/all_sequences.fa
TMPVCF=$OUTVCF.temp
${TABIX} --print-header $INVCF $CHROM | sed 's/ID=AD,Number=./ID=AD,Number=R/' | sed 's/reads with MQ=255 or/reads with MQ equals 255 or/' | ${VT} decompose -s - | ${VT} normalize -r $REF - | ${VT} uniq - | bgzip -c > $TMPVCF
${TABIX} -p vcf -f $TMPVCF && mv $TMPVCF.tbi $OUTVCF.tbi && mv $TMPVCF $OUTVCF
