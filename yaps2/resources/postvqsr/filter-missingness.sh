#!/bin/bash

BIO_1662=/gscmnt/gc2802/halllab/idas/jira/BIO-1662
SCRIPT=${BIO_1662}/bin/identify-missingness
TABIX=/gsc/bin/tabix
GZIP_CMD=/bin/gzip

DBSNP=/gscuser/kmeltzst/gscmnt/reference_files/gotcloud.ref/hapmap_3.3.b37.sites.vcf.gz

INVCF=$1
OUTVCF=$2
STATS=$3

set -o xtrace
${SCRIPT} --stats=${STATS} --db=${DBSNP} --missing-threshold=2.0 ${INVCF} | bgzip -c > ${OUTVCF} \
    && ${TABIX} -p vcf -f ${OUTVCF} \
    && ${GZIP_CMD} ${STATS};
set +o xtrace

