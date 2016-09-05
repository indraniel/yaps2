#!/bin/bash

# modified from @aregier
#  ~aregier/scratch/dlarson-gatk-scripts/run_decompose.sh

BIO_1662=/gscmnt/gc2802/halllab/idas/jira/BIO-1662
VT=${BIO_1662}/vendor/local/bin/vt-0.5

INVCF=$1
OUTVCF=$2

if [ -a $OUTVCF ]
then
    exit 0
fi

REF=/gscmnt/gc2719/halllab/genomes/human/GRCh37/1kg_phase1/human_g1k_v37.fasta
TMPVCF=$OUTVCF.temp
zcat $INVCF | sed 's/ID=AD,Number=./ID=AD,Number=R/' | sed 's/reads with MQ=255 or/reads with MQ equals 255 or/' | ${VT} decompose -s - | ${VT} normalize -r $REF - | ${VT} uniq - | bgzip -c > $TMPVCF
tabix -p vcf -f $TMPVCF && mv $TMPVCF.tbi $OUTVCF.tbi && mv $TMPVCF $OUTVCF
