#!/bin/bash

PLINK=/gscuser/dlarson/bin/plink

IN=$1
EXTRACT=$2
OUT=$3

plink_cmd="${PLINK} --bfile ${IN} --extract ${EXTRACT} --make-bed --out ${OUT}"

echo ${plink_cmd}
eval ${plink_cmd}
