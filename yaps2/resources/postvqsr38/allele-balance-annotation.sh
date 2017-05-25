#!/bin/bash

set -ueo pipefail

BGZIP=/gscmnt/gc2802/halllab/idas/software/local/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/local/bin/tabix
BCFTOOLS=/gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4

PYTHON=$1
SCRIPT=$2
INVCF=$3
OUTVCF=$4
CHROM=$5

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

function is_empty_vcf {
    local invcf=$1

    local count=$(${BCFTOOLS} view -H ${invcf} | head -n 1000 | wc -l)
    if [[ "${count}" -gt "0" ]]; then
        return 1
    else
        return 0
    fi
}

function copy_over_vcf {
    local invcf=$1
    local outvcf=$2

    cp -v ${invcf} ${outvcf}
    cp -v ${invcf}.tbi ${outvcf}.tbi
}

function allele_balance_annotation {
    if [ -e $OUTVCF ]
    then
        exit 0
    fi

    TMPVCF=$OUTVCF.temp

    set -o xtrace
    ${PYTHON} ${SCRIPT} \
        --region="${CHROM}" \
        ${INVCF} \
        | ${BGZIP} -c > ${TMPVCF} \
        && ${TABIX} -p vcf -f ${TMPVCF} \
        && mv $TMPVCF.tbi $OUTVCF.tbi \
        && mv $TMPVCF $OUTVCF ;
    set +o xtrace
}

function main {
    if is_empty_vcf ${INVCF}; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${INVCF} ${OUTVCF} ;
    else
        log "Performing allele balance annotations"
        allele_balance_annotation ;
    fi
}

main ;
