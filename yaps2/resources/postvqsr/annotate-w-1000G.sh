#!/bin/bash

BCFTOOLS=/gsc/bin/bcftools1.2
TABIX=/gsc/bin/tabix

BIO_1984=/gscmnt/gc2802/halllab/idas/jira/BIO-1984

# the 1000G annotation file
KGVCF=${BIO_1984}/data/manual/create-1000G-reformatted-af-annotations/
KGVCF+=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.decompose.normalize.reheader.w_ids.reformatted_pop_af.vcf.gz

INVCF=$1
OUTVCF=$2

${BCFTOOLS} annotate \
    -a ${KGVCF} \
    -c ID,INFO/1KG_EAS_AF,INFO/1KG_EUR_AF,INFO/1KG_AFR_AF,INFO/1KG_AMR_AF,INFO/1KG_SAS_AF \
    -O z \
    -o ${OUTVCF} \
    ${INVCF} \
    && ${TABIX} -p vcf -f ${OUTVCF}
