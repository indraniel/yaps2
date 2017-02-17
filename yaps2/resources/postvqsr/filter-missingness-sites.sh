#!/bin/bash

TABIX=/gsc/bin/tabix
GZIP_CMD=/bin/gzip
BGZIP=/gsc/bin/bgzip1.2.1

DBSNP=/gscuser/kmeltzst/gscmnt/reference_files/gotcloud.ref/hapmap_3.3.b37.sites.vcf.gz

PYTHON=$1
SCRIPT=$2
INVCF=$3
OUTVCF=$4
STATS=$5

set -o xtrace
${PYTHON} ${SCRIPT} \
    --soft \
    --stats=${STATS} \
    --db=${DBSNP} \
    --missing-threshold=2.0 \
    ${INVCF} \
    | ${BGZIP} -c > ${OUTVCF} \
    && ${TABIX} -p vcf -f ${OUTVCF} \
    && ${GZIP_CMD} ${STATS};
set +o xtrace
