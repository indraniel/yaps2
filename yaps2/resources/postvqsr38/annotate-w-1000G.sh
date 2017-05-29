#!/bin/bash

set -eo pipefail

# http://stackoverflow.com/questions/9893667/is-there-a-way-to-write-a-bash-function-which-aborts-the-whole-execution-no-mat
trap "exit 1" TERM
export TOP_PID=$$

JAVA=/gapp/x64linux/opt/java/jdk/jdk1.8.0_60/bin/java
PICARD=/gscmnt/gc2802/halllab/idas/software/picard/picard.2.9.0.jar

PYTHON=$(which python) # if run inside yaps2 pipeline, then should be getting the virtualenv python

BGZIP=/gscmnt/gc2802/halllab/idas/software/local/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/local/bin/tabix
BCFTOOLS=/gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4

function die {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "[ ${timestamp} ] ERROR: $@" >&2
    kill -s TERM ${TOP_PID}
}

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

function run_cmd {
    local cmd=$1
    log "EXEC: ${cmd}"
    eval "${cmd}"
    if [[ $? -ne 0 ]]; then
        die "[err] Problem running command: ${cmd} !"
        exit 1;
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

function vcf_subtract_samples {
    local vcf=$1

    log "(vcf_subtract_samples) calling bcftools subprocess"
    local cmd="
    cat <(${BCFTOOLS} view -h ${vcf} | head -n -1) \
        <(${BCFTOOLS} view -h ${vcf} | tail -n 1 | cut -f1-8) \
        <(${BCFTOOLS} view -H ${vcf} | cut -f1-8)
    "
    run_cmd "${cmd}"
}

function vcf_add_samples {
    local no_samples_vcf=$1
    local samples_vcf=$2

    log "(vcf_add_samples) calling bcftools subprocess"
    local cmd="
    cat <(${BCFTOOLS} view -h ${no_samples_vcf} | head -n -1) \
        <(${BCFTOOLS} view -h ${samples_vcf} | tail -n 1 ) \
        <(paste <(${BCFTOOLS} view -H ${no_samples_vcf}) \
                <(${BCFTOOLS} view -H ${samples_vcf} | cut -f9-))
    "
    run_cmd "${cmd}"
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

function prune_samples_on_b38_vcf {
    local invcf=$1
    local outdir=${2:-$(dirname ${invcf})}
    local outvcf=${outdir}/b38.nosamples.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting prune_samples_on_b38_vcf"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

	local cmd1="vcf_subtract_samples ${invcf} | ${BGZIP} -c > ${tmpvcf}"
	run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}
    echo ${outvcf}
}

function run_liftover_hg19 {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/hg19.vcf.gz
    local reject=${outdir}/hg19.unmapped.vcf.gz
    local logfile=${outdir}/picard.hg19.log

    if [[ -e "${outvcf}" ]] && [[ -e "${outvcf}.tbi" ]] && [[ -e "${logfile}" ]] \
        && grep -q 'picard.vcf.LiftoverVcf done' ${logfile} ; then
        log "shortcutting run_liftover_hg19"
        echo ${outvcf}
        return 0;
    fi

    local chain=/gscmnt/gc2802/halllab/aregier/jira/BIO-2228/hg38ToHg19.over.chain.gz
    local reference=/gscmnt/gc2719/halllab/genomes/human/GRCh37/hg19_ucsc/hg19.fa

    local cmd1="
    ${JAVA} \
        -Xmx16g \
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

    if [[ -e "${outvcf}" ]] && [[ -e "${outvcf}.tbi" ]] && [[ -e "${logfile}" ]] \
        && grep -q 'picard.vcf.LiftoverVcf done' ${logfile} ; then
        log "shortcutting run_liftover_grc37"
        echo ${outvcf}
        return 0;
    fi

    local chain=/gscmnt/gc2802/halllab/aregier/jira/BIO-2228/hg19ToGRCh37.over.chain.gz
    local reference=/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa

    local cmd1="
    ${JAVA} \
        -Xmx16g \
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

function run_1kg_annotation {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b37-1kg-annotation.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_1kg_annotation"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # the 1000G annotation file
    local BIO_1984=/gscmnt/gc2802/halllab/idas/jira/BIO-1984
    local KGVCF=${BIO_1984}/data/manual/create-1000G-reformatted-af-annotations/
    KGVCF+=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.decompose.normalize.reheader.w_ids.reformatted_pop_af.vcf.gz

    local cmd1="
    ${BCFTOOLS} annotate \
        -a ${KGVCF} \
        -c INFO/1KG_EAS_AF,INFO/1KG_EUR_AF,INFO/1KG_AFR_AF,INFO/1KG_AMR_AF,INFO/1KG_SAS_AF \
        -O z \
        -o ${tmpvcf} \
        ${invcf} \
    "
    run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}

    echo ${outvcf}
}

function integrate_b37_annotations_to_b38 {
    local integrate_script=$1
    local b37_vcf=$2
    local b38_vcf=$3
    local annotation_type=$4

    local outdir=$(dirname ${b38_vcf})
    local outvcf=${outdir}/b38.1kg.nosamples.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting integrate_b37_annotations_to_b38"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    local cmd1="
    ${PYTHON} ${integrate_script} \
        --b38-vcf=${b38_vcf}  \
        --annotated-b37-vcf=${b37_vcf} \
        --auto-fill \
        --annotation-type=${annotation_type} \
    | ${BGZIP} -c \
    > ${tmpvcf}
    "
    run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}
    echo ${outvcf}
}

function add_samples_on_b38_annotated_vcf {
	local original_b38_input_vcf=$1
    local b38_annotated_no_samples_vcf=$2
    local final_b38_output_vcf=$3

    if [[ -e "${final_b38_output_vcf}" ]]; then
        log "shortcutting add_samples_on_b38_cadd_vcf"
        return 0;
    fi

    local tmpvcf=${final_b38_output_vcf}.tmp

    local cmd1="
	vcf_add_samples ${b38_annotated_no_samples_vcf} ${original_b38_input_vcf} \
        | ${BGZIP} -c \
        > ${tmpvcf}
    "
    run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${final_b38_output_vcf}
}

function annotate_vcf {
    local b38_invcf=$1
    local b38_outvcf=$2
    local integrate_script=$3

    local scratch_dir=$(dirname ${b38_outvcf})/scratch
    mkdir -p ${scratch_dir}

    log "Remove samples on b38 input vcf"
    local b38_invcf_no_samples=$(prune_samples_on_b38_vcf ${b38_invcf} ${scratch_dir})
    log "Entering liftOver hg19"
    local hg19_vcf=$(run_liftover_hg19 ${b38_invcf_no_samples})
    log "Entering liftOver GRCh37"
    local grc37_vcf=$(run_liftover_grc37 ${hg19_vcf})
    log "Entering run_1kg_annotation"
    local b37_1kg_vcf=$(run_1kg_annotation ${grc37_vcf})
    log "Entering integrate b37 cadd annotations back to b38"
    local b38_anno_vcf=$(integrate_b37_annotations_to_b38 ${integrate_script} ${b37_1kg_vcf} ${b38_invcf_no_samples} '1000G')
    log "Add samples on b38 cadd annotated vcf"
    add_samples_on_b38_annotated_vcf ${b38_invcf} ${b38_anno_vcf} ${b38_outvcf}
}

function main {
    local invcf=$1
    local outvcf=$2
    local integrate_script=$3

    if is_empty_vcf ${invcf} ; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${invcf} ${outvcf} ;
    else
        annotate_vcf ${invcf} ${outvcf} ${integrate_script};
    fi

    log 'All Done'
}

INVCF=$1
OUTVCF=$2
INTEGRATE_SCRIPT=$3

main ${INVCF} ${OUTVCF} ${INTEGRATE_SCRIPT} ;
