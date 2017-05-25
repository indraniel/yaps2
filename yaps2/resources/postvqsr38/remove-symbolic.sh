#!/bin/bash

set -ueo pipefail

BCFTOOLS=/gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4
TABIX=/gscmnt/gc2802/halllab/idas/software/local/bin/tabix

INVCF=$1
OUTVCF=$2

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

function remove_symbolic {
    if [ -e $OUTVCF ]
    then
        exit 0
    fi

    TMPVCF=$OUTVCF.temp
    ${BCFTOOLS} view -e '%TYPE="other" || ALT="*"' $INVCF --output-type z --output-file $TMPVCF \
        && ${TABIX} -p vcf -f $TMPVCF && mv $TMPVCF.tbi $OUTVCF.tbi && mv $TMPVCF $OUTVCF
}

function main {
    if is_empty_vcf ${INVCF}; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${INVCF} ${OUTVCF} ;
    else
        remove_symbolic ;
    fi
}

main ;
