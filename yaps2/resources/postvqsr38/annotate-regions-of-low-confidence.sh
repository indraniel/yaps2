#!/bin/bash

set -eo pipefail

# http://stackoverflow.com/questions/9893667/is-there-a-way-to-write-a-bash-function-which-aborts-the-whole-execution-no-mat
trap "exit 1" TERM
export TOP_PID=$$

BGZIP=/gscmnt/gc2802/halllab/idas/software/local/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/local/bin/tabix
BCFTOOLS=/gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4

function join_by { local IFS="$1"; shift; echo "$*"; }

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

function run_LCR_annotation {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b38-LCR-annotation.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_LCR_annotation"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # the LCR annotation file
    local base=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations
    local anno_vcf=${base}/LCR-hs38.bed.gz

    if [[ ! -e ${anno_vcf} ]]; then
        die "Did not find annotation vcf: '${anno_vcf}'"
    fi

    # annotation columns of interest
    local -a anno_cols=(
        'CHROM'
        'FROM'
        'TO'
    )
    local annotations=$(join_by , "${anno_cols[@]}")
    local hdr_msg="Indicates that the variant is in a Low Confidence Region (LCR)"

    local cmd1="
    ${BCFTOOLS} annotate \
        -a ${anno_vcf} \
        -c ${annotations} \
        -h <(echo '##INFO=<ID=LCR,Number=0,Type=Flag,Description="${hdr_msg}">') \
        -m LCR \
        -O z \
        -o ${tmpvcf} \
        ${invcf} \
    "
    run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}

    echo ${outvcf}
}

function run_centromeres_annotation {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b38-centromeres-annotation.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_centromeres_annotation"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # the centromeres annotation file
    local base=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations
    local anno_vcf=${base}/centromeres.bed.gz

    if [[ ! -e ${anno_vcf} ]]; then
        die "Did not find annotation vcf: '${anno_vcf}'"
    fi

    # annotation columns of interest
    local -a anno_cols=(
        'CHROM'
        'FROM'
        'TO'
    )
    local annotations=$(join_by , "${anno_cols[@]}")
    local hdr_msg="Variant is in a centromere region (CEN)"

    local cmd1="
    ${BCFTOOLS} annotate \
        -a ${anno_vcf} \
        -c ${annotations} \
        -h <(echo '##INFO=<ID=CEN,Number=0,Type=Flag,Description="${hdr_msg}">') \
        -m CEN \
        -O z \
        -o ${tmpvcf} \
        ${invcf} \
    "
    run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}

    echo ${outvcf}
}

function run_segdups_annotation {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b38-segdups-annotation.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_segdups_annotation"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # the segdups annotation file
    local base=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations
    local anno_vcf=${base}/segdups.bed.gz

    if [[ ! -e ${anno_vcf} ]]; then
        die "Did not find annotation vcf: '${anno_vcf}'"
    fi

    # annotation columns of interest
    local -a anno_cols=(
        'CHROM'
        'FROM'
        'TO'
        'SEGDUPS'
    )
    local annotations=$(join_by , "${anno_cols[@]}")
    local hdr_msg="segdups region"

    local cmd1="
    ${BCFTOOLS} annotate \
        -a ${anno_vcf} \
        -c ${annotations} \
        -h <(echo '##INFO=<ID=SEGDUPS,Number=1,Type=Integer,Description="${hdr_msg}">') \
        -O z \
        -o ${tmpvcf} \
        ${invcf} \
    "
    run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}

    echo ${outvcf}
}

function run_satellite_annotation {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b38-satellite-annotation.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_satellite_annotation"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # the satellite annotation file
    local base=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations
    local anno_vcf=${base}/satellite.hg38.bed.gz

    if [[ ! -e ${anno_vcf} ]]; then
        die "Did not find annotation vcf: '${anno_vcf}'"
    fi

    # annotation columns of interest
    local -a anno_cols=(
        'CHROM'
        'FROM'
        'TO'
        'SATELLITE'
    )
    local annotations=$(join_by , "${anno_cols[@]}")
    local hdr_msg="satellite region"

    local cmd1="
    ${BCFTOOLS} annotate \
        -a ${anno_vcf} \
        -c ${annotations} \
        -h <(echo '##INFO=<ID=SATELLITE,Number=1,Type=String,Description="${hdr_msg}">') \
        -O z \
        -o ${tmpvcf} \
        ${invcf} \
    "
    run_cmd "${cmd1}"
	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}

    echo ${outvcf}
}

function run_full_pipeline_annotation {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b38-full-lcr-annotation.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_full_pipeline_annotation"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # the annotation files
    local base=/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/annotations
    local lcr_vcf=${base}/LCR-hs38.bed.gz
    local centromeres_vcf=${base}/centromeres.bed.gz
    local segdups_vcf=${base}/segdups.bed.gz
    local satellites_vcf=${base}/satellite.hg38.bed.gz

    for vcf in ${lcr_vcf} ${centromeres_vcf} ${segdups_vcf} ${satellites_vcf}; do
        if [[ ! -e ${vcf} ]]; then
            die "Did not find annotation vcf: '${vcf}'"
        fi
    done

    # setup LCR annotation params
    local lcr_label="LCR"
    local lcr_number="0"
    local lcr_type="Flag"
    local -a anno_lcr_cols=( 'CHROM' 'FROM' 'TO' )
    local lcr_annotations=$(join_by , "${anno_lcr_cols[@]}")
    local lcr_hdr_msg="Indicates that the variant is in a Low Confidence Region (LCR)"
    local lcr_header="'##INFO=<ID=${lcr_label},Number=${lcr_number},Type=${lcr_type},Description=\"${lcr_hdr_msg}\">'"

    # setup centromeres annotation params
    local centromeres_label="CEN"
    local centromeres_number="0"
    local centromeres_type="Flag"
    local centromeres_annotations=$(join_by , "${anno_lcr_cols[@]}")
    local centromeres_hdr_msg="Variant is in a centromere region (CEN)"
    local centromeres_header="'##INFO=<ID=${centromeres_label},Number=${centromeres_number},Type=${centromeres_type},Description=\"${centromeres_hdr_msg}\">'"

    # setup segdups annotation params
    local segdups_label="SEGDUPS"
    local segdups_number="1"
    local segdups_type="Integer"
    local -a anno_segdups_cols=( 'CHROM' 'FROM' 'TO' 'SEGDUPS')
    local segdups_annotations=$(join_by , "${anno_segdups_cols[@]}")
    local segdups_hdr_msg="segdups region"
    local segdups_header="'##INFO=<ID=${segdups_label},Number=${segdups_number},Type=${segdups_type},Description=\"${segdups_hdr_msg}\">'"

    # setup satellites annotation params
    local satellites_label="SATELLITE"
    local satellites_number="1"
    local satellites_type="String"
    local -a anno_satellites_cols=( 'CHROM' 'FROM' 'TO' 'SATELLITE')
    local satellites_annotations=$(join_by , "${anno_satellites_cols[@]}")
    local satellites_hdr_msg="satellite region"
    local satellites_header="'##INFO=<ID=${satellites_label},Number=${satellites_number},Type=${satellites_type},Description=\"${satellites_hdr_msg}\">'"

    local cmd1="
    ${BCFTOOLS} annotate \
        -a ${lcr_vcf} \
        -c ${lcr_annotations} \
        -h <(echo "${lcr_header}") \
        -m ${lcr_label} \
        -O u \
        ${invcf} \
    | ${BCFTOOLS} annotate \
        -a ${centromeres_vcf} \
        -c ${centromeres_annotations} \
        -h <(echo "${centromeres_header}") \
        -m ${centromeres_label} \
        -O u \
        - \
    | ${BCFTOOLS} annotate \
        -a ${segdups_vcf} \
        -c ${segdups_annotations} \
        -h <(echo "${segdups_header}") \
        -O u \
        - \
    | ${BCFTOOLS} annotate \
        -a ${satellites_vcf} \
        -c ${satellites_annotations} \
        -h <(echo "${satellites_header}") \
        -O z \
        -o ${tmpvcf} \
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

    local scratch_dir=$(dirname ${b38_outvcf})/scratch
    mkdir -p ${scratch_dir}

    log "Remove samples on b38 input vcf"
    local b38_invcf_no_samples=$(prune_samples_on_b38_vcf ${b38_invcf} ${scratch_dir})

#    log "Entering run_LCR_annotation"
#    local LCR_anno_vcf=$(run_LCR_annotation ${b38_invcf_no_samples})
#    log "Entering run_centromeres_annotation"
#    local centromeres_anno_vcf=$(run_centromeres_annotation ${LCR_anno_vcf})
#    log "Entering run_segdups_annotation"
#    local segdups_anno_vcf=$(run_segdups_annotation ${LCR_anno_vcf})
#    log "Entering run_satellite_annotation"
#    local satellite_anno_vcf=$(run_satellite_annotation ${LCR_anno_vcf})
#    log "Add samples on b38 cadd annotated vcf"
#    # the final annotated vcf file in the pipeline
#    local b38_anno_vcf=${satellite_anno_vcf} 

    log "Entering run_full_pipeline_annotation"
    local b38_anno_vcf=$(run_full_pipeline_annotation ${b38_invcf_no_samples})

    add_samples_on_b38_annotated_vcf ${b38_invcf} ${b38_anno_vcf} ${b38_outvcf}
}

function main {
    local invcf=$1
    local outvcf=$2

    if is_empty_vcf ${invcf} ; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${invcf} ${outvcf} ;
    else
        annotate_vcf ${invcf} ${outvcf} ;
    fi

    log 'All Done'
}

INVCF=$1
OUTVCF=$2

main ${INVCF} ${OUTVCF} ;
