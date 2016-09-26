#!/bin/bash

set -ueo pipefail

BCFTOOLS=/gsc/bin/bcftools1.2

FAM_FILE=$1
VCF=$2
FILTER=${3-}

while read FAM IND DAD MOM stuff
do
    if [[ "$DAD" != "0" && "$MOM" != "0" ]]; then
        OPT=""
        if [[ "$FILTER" ]]; then
            OPT="-f .,PASS"
        fi

        RESULT=$(${BCFTOOLS} view -s $IND,$DAD,$MOM $OPT $VCF | ${BCFTOOLS} view --min-ac 1 - | ${BCFTOOLS} view -s $DAD,$MOM - | ${BCFTOOLS} view -g hom - | grep -v 'SECONDARY' | { grep -v '^#' || true; } | wc -l)
        echo -e "$FAM\t$RESULT"
    fi
done < $FAM_FILE
