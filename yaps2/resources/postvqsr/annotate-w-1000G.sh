#!/bin/bash

BCFTOOLS=/gsc/bin/bcftools1.2
TABIX=/gsc/bin/tabix

BIO_1662=/gscmnt/gc2802/halllab/idas/jira/BIO-1662

# the 1000G annotation file
KGVCF=${BIO_1662}/data/derived/FinnMetSeq-WGS/1000-regenerate-1000G-annotation-vcf/
KGVCF+=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.decompose.normalize.reheader.w_ids.vcf.gz

INVCF=$1
OUTVCF=$2

${BCFTOOLS} annotate -a ${KGVCF} -c ID -O z -o ${OUTVCF} ${INVCF} \
    && ${TABIX} -p vcf -f ${OUTVCF}
