#!/bin/bash

BGZIP=/gscmnt/gc2802/halllab/idas/software/local/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/local/bin/tabix

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
