#!/bin/bash

BCFTOOLS=/gsc/bin/bcftools1.2
TABIX=/gsc/bin/tabix

BIO_1984=/gscmnt/gc2802/halllab/idas/jira/BIO-1984

# the ExAC annotation file
KGVCF=${BIO_1984}/data/manual/create-ExAC-reformatted-af-annotations/
KGVCF+=ExAC.r0.3.1.sites.simpleAF.multiallelic.reformatted_pop_af.vcf.gz

INVCF=$1
OUTVCF=$2

${BCFTOOLS} annotate \
    -a ${KGVCF} \
    -c INFO/EXAC_AFR_AF,INFO/EXAC_AMR_AF,INFO/EXAC_Adj_AF,INFO/EXAC_EAS_AF,INFO/EXAC_FIN_AF,INFO/EXAC_NFE_AF,INFO/EXAC_OTH_AF,INFO/EXAC_SAS_AF \
    -O z \
    -o ${OUTVCF} \
    ${INVCF} \
    && ${TABIX} -p vcf -f ${OUTVCF}
