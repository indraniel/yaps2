#!/bin/bash

PLINK=/gscuser/dlarson/bin/plink
PERL=/usr/bin/perl

IN=$1
MERGE_LIST=$2
OUT=$3

plink_cmd="${PLINK} --bfile ${IN} --merge-list ${MERGE_LIST} --recode tab --out ${OUT}"

# address downstream eigenstrat issue: 'convertf decides to “ignore” all my samples. Why?'
# see: https://www.hsph.harvard.edu/alkes-price/eigensoft-frequently-asked-questions/
perl_cmd="${PERL} -p -i -e '\$_ =~ s/\t-9\t/\t1\t/' ${OUT}.ped"

full_cmd="${plink_cmd} && ${perl_cmd}"

echo ${full_cmd}
eval ${full_cmd}
