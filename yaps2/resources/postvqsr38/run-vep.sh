#!/bin/bash

set -eo pipefail

# this script is meant to be run inside docker image willmclaren/ensembl-vep:release_88
# see: https://hub.docker.com/r/willmclaren/ensembl-vep/
#      https://github.com/Ensembl/ensembl-vep/blob/release/88/docker/Dockerfile

function join_by { local IFS="$1"; shift; echo "$*"; }

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

function run_cmd {
    local cmd=$1
    log "EXEC: ${cmd}"
    eval "${cmd}"
    if [[ $? -ne 0 ]]; then
        log "[err] Problem running command: ${cmd} !"
        exit 1
    fi
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

function tabix_and_finalize_vcf {
	local tmpvcf=$1
	local finalvcf=$2
    local cmd="
    ${TABIX} -p vcf -f ${tmpvcf} \
        && mv ${tmpvcf}.tbi ${finalvcf}.tbi \
        && mv ${tmpvcf} ${finalvcf}
    "
    run_cmd "${cmd}"
}

function run_vep {
    local invcf=$1
    local outdir=$2
    local outvcf=${outdir}/vep.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_vep"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=$outvcf.tmp

    # see BIO-2229 for vep parameter details
    local vep=/home/vep/src/ensembl-vep/vep
    local vep_cache=/gscmnt/gc2719/halllab/genomes/human/GRCh38/.vep
    local vep_fields="Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE"

    # ensure that the correct perl and libraries are being used
    local -a extra_perl_libs=(
        '/home/vep/src/ensembl-vep'
        '/home/vep/src/bioperl-live-release-1-6-924'
        '/usr/local/lib/x86_64-linux-gnu/perl/5.22.1'
        '/usr/local/share/perl/5.22.1'
        '/usr/share/perl/5.22.1'
    )
    local perl_libs=$(join_by : "${extra_perl_libs[@]}")

    local perl_path=/gscmnt/gc2802/halllab/idas/software/std-perl
    local htslib_path=/home/vep/src/htslib
    export PERL5LIB=${perl_libs}
    export PATH=${perl_path}:${htslib_path}:${PATH}

    local cmd1="
    zcat ${invcf} \
        | ${vep} \
            --force_overwrite \
            --offline \
            --fork 12 \
            --cache \
            --dir_cache ${vep_cache} \
            --sift b  \
            --polyphen b  \
            --species homo_sapiens  \
            --symbol   \
            --numbers \
            --biotype \
            --total_length \
            -o STDOUT  \
            --format vcf \
            --vcf \
            --fields ${vep_fields} \
            --no_stats \
        | ${BGZIP} -c \
        > $tmpvcf ;
    "
    run_cmd "${cmd1}"

	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}

    echo ${outvcf}
}

function main {
    local invcf=$1
    local outvcf=$2

    if is_empty_vcf ${invcf} ; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${invcf} ${outvcf} ;
    else
        local outdir=$(dirname ${outvcf})
        log "Entering run_vep"
        local vepvcf=$(run_vep ${invcf} ${outdir})
        local cmd2="
            mv -v ${vepvcf}.tbi ${outvcf}.tbi \
            && mv -v ${vepvcf} ${outvcf}
        "
        run_cmd "${cmd2}"
    fi

    log 'All Done'
}

INVCF=$1
OUTVCF=$2

main ${INVCF} ${OUTVCF};
