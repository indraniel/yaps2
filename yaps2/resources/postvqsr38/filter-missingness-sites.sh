#!/bin/bash

TABIX=/gsc/bin/tabix
GZIP_CMD=/bin/gzip
BGZIP=/gsc/bin/bgzip1.2.1

PYTHON=$1
SCRIPT=$2
INVCF=$3
OUTVCF=$4
CHROM=$5

set -o xtrace
${PYTHON} ${SCRIPT} \
    --soft \
    --missing-threshold=2.0 \
    --region="${CHROM}" \
    ${INVCF} \
    | ${BGZIP} -c > ${OUTVCF} \
    && ${TABIX} -p vcf -f ${OUTVCF} ;
set +o xtrace
