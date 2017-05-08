#!/bin/bash

set -eo pipefail

# this script is meant to be run inside docker image willmclaren/ensembl-vep:release_88
# see: https://hub.docker.com/r/willmclaren/ensembl-vep/
#      https://github.com/Ensembl/ensembl-vep/blob/release/88/docker/Dockerfile

AWK=/usr/bin/awk

JAVA=/gapp/x64linux/opt/java/jdk/jdk1.8.0_60/bin/java
PICARD=/gscmnt/gc2802/halllab/idas/software/picard/picard.2.9.0.jar

PYTHON=$(which python) # if run inside yaps2 pipeline, then should be getting the virtualenv python

BGZIP=/gscmnt/gc2802/halllab/idas/software/vep/local/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/vep/local/bin/tabix
BCFTOOLS=/gscmnt/gc2802/halllab/idas/software/vep/local/bcftools1.4

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

function subtract_vcf_samples {
    local invcf=$1
    local new_vcf_name=$2
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/${new_vcf_name}
    local tmpvcf=$outvcf.tmp

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting subtract_vcf_samples (${new_vcf_name})"
        echo ${outvcf}
        return 0;
    fi

    # add the header to the new test vcf (except for the last header line)
    cmd1="${BCFTOOLS} view -h ${invcf} | head -n -1 > ${tmpvcf}"
    run_cmd "${cmd1}"

    # correctly add the last header line
    cmd2="${BCFTOOLS} view -h ${invcf} | tail -n 1 | cut -f 1-8 >>${tmpvcf}"
    run_cmd "${cmd2}"

    # add the first few samples / variants
    cmd3="${BCFTOOLS} view -H ${invcf} | cut -f 1-8 >>${tmpvcf}"
    run_cmd "${cmd3}"

    # tabix & bgzip the file
    local cmd4="
    ${TABIX} -p vcf -f ${tmpvcf} \
        && mv ${tmpvcf}.tbi ${outvcf}.tbi \
        && mv ${tmpvcf} ${outvcf}
    "
    run_cmd "${cmd4}"
    echo ${outvcf}
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

    local cmd2="
    ${TABIX} -p vcf -f ${tmpvcf} \
        && mv ${tmpvcf}.tbi ${outvcf}.tbi \
        && mv ${tmpvcf} ${outvcf}
    "
    run_cmd "${cmd2}"

    echo ${outvcf}
}

function run_liftover_hg19 {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/hg19.vcf.gz
    local reject=${outdir}/hg19.unmapped.vcf.gz
    local logfile=${outdir}/picard.hg19.log

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_liftover_hg19"
        echo ${outvcf}
        return 0;
    fi

    local chain=/gscmnt/gc2802/halllab/aregier/jira/BIO-2228/hg38ToHg19.over.chain.gz
    local reference=/gscmnt/gc2719/halllab/genomes/human/GRCh37/hg19_ucsc/hg19.fa

    local cmd1="
    ${JAVA} \
        -Xmx8g \
        -jar ${PICARD} \
        LiftoverVcf \
        I=${invcf} \
        O=${outvcf} \
        C=${chain} \
        REJECT=${reject} \
        R=${reference} \
        WRITE_ORIGINAL_POSITION=true \
        >${logfile} 2>&1
    "

    run_cmd "${cmd1}"
    echo ${outvcf}
}

function run_liftover_grc37 {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/grc37.vcf.gz
    local reject=${outdir}/grc37.unmapped.vcf.gz
    local logfile=${outdir}/picard.grc37.log

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_liftover_grc37"
        echo ${outvcf}
        return 0;
    fi

    local chain=/gscmnt/gc2802/halllab/aregier/jira/BIO-2228/hg19ToGRCh37.over.chain.gz
    local reference=/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa

    local cmd1="
    ${JAVA} \
        -Xmx8g \
        -jar ${PICARD} \
        LiftoverVcf \
        I=${invcf} \
        O=${outvcf} \
        C=${chain} \
        REJECT=${reject} \
        R=${reference} \
        WRITE_ORIGINAL_POSITION=false \
        >${logfile} 2>&1
    "

    run_cmd "${cmd1}"
    echo ${outvcf}
}

function run_cadd {
    local invcf=$1
    local outdir=$(dirname ${invcf})

    local tsv=${outdir}/cadd-annotation.tsv.gz

    if [[ -e "${tsv}" ]]; then
        log "shortcutting run_cadd"
        echo ${tsv}
        return 0;
    fi

    local tmptsv=${tsv}.tmp

    # ensure that the correct programs and libraries are being used
    local -a extra_perl_libs=(
        '/home/vep/src/ensembl-vep'
        '/home/vep/src/bioperl-live-release-1-6-924'
        '/usr/local/lib/x86_64-linux-gnu/perl/5.22.1'
        '/usr/local/share/perl/5.22.1'
        '/usr/share/perl/5.22.1'
    )
    local perl_libs=$(join_by : "${extra_perl_libs[@]}")
    local vep_path=/home/vep/src/ensembl-vep/vep
    local htslib_path=/home/vep/src/htslib
    local perl_path=/gscmnt/gc2802/halllab/idas/software/std-perl
    local awk_path=/gscmnt/gc2802/halllab/idas/software/std-awk
    export PERL5LIB=${perl_libs}
    export PATH=${vep_path}:${perl_path}:${awk_path}:${htslib_path}:${PATH}

    local cadd_dir=/gscmnt/gc2802/halllab/idas/software/cadd-1.2/CADD_v1.2
    local cadd_score_script=${cadd_dir}/bin/score.sh

    local cmd="/bin/bash ${cadd_score_script} ${invcf} ${tmptsv} && mv ${tmptsv} ${tsv}"
    run_cmd "${cmd}"

    echo ${tsv}
}

function paste_cadd {
    local tsvfile=$1
    local vepvcf=$2
    local origvcf=$3
    local merge_cadd_script=$4

    local outdir=$(dirname ${vepvcf})
    local outvcf=${outdir}/grc37.vep.cadd.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting paste_cadd"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # core of "paste_cadd.sh"
    local cmd1="
    cat <(${PYTHON} ${merge_cadd_script} -H <(zcat ${vepvcf}) <(zcat ${tsvfile})) \
        <(paste <(${PYTHON} ${merge_cadd_script} --no-header <(zcat ${vepvcf}) <(zcat ${tsvfile})) \
                <(zcat ${origvcf} | ${AWK} '{if(\$5!=\".\") print \$0}' | grep -v '^##' | cut -f9- )) \
        | ${BGZIP} -c \
        > ${tmpvcf}
    "
    run_cmd "${cmd1}"

    local cmd2="
    ${TABIX} -p vcf -f ${tmpvcf} \
        && mv ${tmpvcf}.tbi ${outvcf}.tbi \
        && mv ${tmpvcf} ${outvcf}
    "
    run_cmd "${cmd2}"

    echo ${outvcf}
}

function reset_to_b38 {
    local invcf=$1
    local outvcf=$2
    local reset_script=$3

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting reset_to_b38"
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    log "Please implement me!"
    exit 1

    echo ${outvcf}
}

function annotate_vcf {
    local invcf=$1
    local outvcf=$2
    local merge_script=$3
    local reset_liftover_script=$4


    local outdir=$(dirname ${outvcf})

    log "Entering run_vep"
    local vepvcf=$(run_vep ${invcf} ${outdir})
    log "Entering liftOver hg19"
    local hg19vcf=$(run_liftover_hg19 ${vepvcf})
    log "Entering liftOver GRCh37"
    local grc37vcf=$(run_liftover_grc37 ${hg19vcf})
    log "Remove samples on GRCh37 liftOver vcf"
    local no_sample_grc37_vcf=$(subtract_vcf_samples ${grc37vcf} 'grc37.nosample.vcf.gz')
    log "Entering run_cadd"
    local tsv=$(run_cadd ${no_sample_grc37vcf})
    log "Entering paste_cadd"
    local b37vcf=$(paste_cadd ${tsv} ${no_sample_grc37_vcf} ${grc37vcf} ${merge_script})
    log "Entering annotate back to b38"
    reset_to_b38 ${b37vcf} ${outvcf} ${reset_liftover_script}
}

function main {
    local invcf=$1
    local outvcf=$2
    local merge_script=$3
    local reset_liftover_script=$4

    if is_empty_vcf ${invcf} ; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${invcf} ${outvcf} ;
    else
        annotate_vcf ${invcf} ${outvcf} ${merge_script} ${reset_liftover_script};
    fi

    log 'All Done'
}

INVCF=$1
OUTVCF=$2
MERGE_SCRIPT=$3
RESET_LIFTOVER_SCRIPT=$4

main ${INVCF} ${OUTVCF} ${MERGE_SCRIPT} ${RESET_LIFTOVER_SCRIPT};
