#!/bin/bash

BCFTOOLS=/gsc/bin/bcftools1.2
BGZIP=/gsc/bin/bgzip1.2.1
TABIX=/gsc/bin/tabix1.2.1

INVCF=$1
OUTVCF=$2
MIN_VQSLOD=$3

TMPVCF=${OUTVCF}.temp

bcftools_cmd="${BCFTOOLS} view -O b -i 'INFO/VQSLOD>=${MIN_VQSLOD}' -f '.,PASS' -m 2 -M 2 -v snps ${INVCF} \
    | ${BCFTOOLS} annotate --set-id +'%CHROM-%POS-%REF-%FIRST_ALT' -O z -o ${TMPVCF} "

tabix_cmd="${TABIX} -p vcf -f ${TMPVCF}"
mv_cmd="mv ${TMPVCF}.tbi ${OUTVCF}.tbi && mv ${TMPVCF} ${OUTVCF}"
full_cmd="${bcftools_cmd} && ${tabix_cmd} && ${mv_cmd}"

echo ${full_cmd}
eval ${full_cmd}
