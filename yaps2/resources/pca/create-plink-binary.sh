#!/bin/bash

PLINK=/gscuser/dlarson/bin/plink

INVCF=$1
OUT=$2

plink_cmd="${PLINK} --vcf ${INVCF} --double-id --make-bed --out ${OUT}"

echo ${plink_cmd}
eval ${plink_cmd}
