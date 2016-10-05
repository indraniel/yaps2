#!/bin/bash

PLINK=/gscuser/dlarson/bin/plink

IN=$1
OUT=$2

MAF=0.05
WINDOW=50
STEP=5
R2=0.3
GENO=0

plink_cmd="${PLINK} --bfile ${IN} --maf ${MAF} --indep-pairwise ${WINDOW} ${STEP} ${R2} --geno ${GENO} --out ${OUT}"

echo ${plink_cmd}
eval ${plink_cmd}
