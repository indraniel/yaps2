#!/bin/bash

set -eo pipefail

# http://stackoverflow.com/questions/9893667/is-there-a-way-to-write-a-bash-function-which-aborts-the-whole-execution-no-mat
trap "exit 1" TERM
export TOP_PID=$$

JAVA=/gapp/x64linux/opt/java/jdk/jdk1.8.0_60/bin/java
PICARD=/gscmnt/gc2802/halllab/idas/software/picard/picard.2.9.0.jar

PYTHON=$(which python) # if run inside yaps2 pipeline, then should be getting the virtualenv python
AWK=/usr/bin/awk
LN=/bin/ln
SORT=/gscmnt/gc2802/halllab/idas/software/local/bin/lh3-sort
UNIQ=/usr/bin/uniq

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

function get_chrom {
    local region=$1
    local chrom=

    if echo "${region}" | grep -q ':'; then
        chrom=$(echo ${region} | ${AWK} -F':' '{ print $1 }' | sed 's/chr//g')
        #interval=$(echo ${region} | ${AWK} -F':' '{ print $2 }')
    else
        chrom=${region}
    fi

	if [[ -z "${chrom}" ]]; then
        die "Could not obtain a chrom from region : ${region}"
	fi

    log "identified chrom '${chrom}' on region '${region}'"
    echo ${chrom}
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

function run_gnomAD_exome_annotation {
    local region=$1
    local invcf=$2
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b37-gnomAD-exome-annotation-final.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_gnomAD_exome_annotation"
        echo ${outvcf}
        return 0;
    fi

    # annotation columns of interest
    local -a anno_cols=(
        'ID'
        'INFO/GNOMAD_EXOME_AC_AFR'
        'INFO/GNOMAD_EXOME_AC_AMR'
        'INFO/GNOMAD_EXOME_AC_ASJ'
        'INFO/GNOMAD_EXOME_AC_EAS'
        'INFO/GNOMAD_EXOME_AC_FIN'
        'INFO/GNOMAD_EXOME_AC_NFE'
        'INFO/GNOMAD_EXOME_AC_OTH'
        'INFO/GNOMAD_EXOME_AC_SAS'
        'INFO/GNOMAD_EXOME_AN_AFR'
        'INFO/GNOMAD_EXOME_AN_AMR'
        'INFO/GNOMAD_EXOME_AN_ASJ'
        'INFO/GNOMAD_EXOME_AN_EAS'
        'INFO/GNOMAD_EXOME_AN_FIN'
        'INFO/GNOMAD_EXOME_AN_NFE'
        'INFO/GNOMAD_EXOME_AN_OTH'
        'INFO/GNOMAD_EXOME_AN_SAS'
        'INFO/GNOMAD_EXOME_AC_raw'
        'INFO/GNOMAD_EXOME_AN_raw'
        'INFO/GNOMAD_EXOME_AC_POPMAX'
        'INFO/GNOMAD_EXOME_AN_POPMAX'
        'INFO/GNOMAD_EXOME_AS_RF'
        'INFO/GNOMAD_EXOME_AS_FilterStatus'
        'INFO/GNOMAD_EXOME_AS_RF_POSITIVE_TRAIN'
        'INFO/GNOMAD_EXOME_AS_RF_NEGATIVE_TRAIN'
    )
    local annotations=$(join_by , "${anno_cols[@]}")

    # the gnomAD exome annotation vcfs base location
    local base=/gscmnt/gc2802/halllab/gnomAD/release-170228/processed/post-vqsr-pipeline/exome

    # Figure out which gnomAD chromosomes vcfs are needed
    local -a gnomAD_chroms=($(${TABIX} --list-chroms ${invcf} | ${SORT} -N | ${UNIQ} | grep -v -P '^(GL|MT)'))
    log "[exome] gnomAD chromosomal vcfs to process: ${gnomAD_chroms[@]} ( items: ${#gnomAD_chroms[@]} )"

    local enter_vcf=${invcf}
    local exit_vcf=

    for chr in ${gnomAD_chroms[@]}; do
        log "Annotating with gnomAD exome chromosome: ${chr}"
        exit_vcf=${outdir}/b37-gnomAD-exome-annotation-c${chr}.vcf.gz

        # determine the annotation vcf
        local anno_vcf=${base}/${chr}/
        anno_vcf+=c${chr}.gnomad.exomes.r2.0.1.sites.decompose.normalized.uniq.namespaced.vcf.gz

#        log "enter_vcf : ${enter_vcf}"
#        log "exit_vcf : ${exit_vcf}"

        if [[ ! -e ${anno_vcf} ]]; then
            die "Did not find annotation vcf: '${anno_vcf}'"
        fi

        # perform the chromosomal annotation
        local tmpvcf=${exit_vcf}.tmp
        local cmd1="
        ${BCFTOOLS} annotate \
            -a ${anno_vcf} \
            -c ${annotations} \
            -O z \
            -o ${tmpvcf} \
            ${enter_vcf} \
        "
        run_cmd "${cmd1}"
        tabix_and_finalize_vcf ${tmpvcf} ${exit_vcf}

        # setup for the next iteration
        enter_vcf=${exit_vcf}
    done

    log "final exome exit_vcf is: ${exit_vcf}"
    local cmd2="mv ${exit_vcf}.tbi ${outvcf}.tbi && mv ${exit_vcf} ${outvcf}"
    run_cmd "${cmd2}"

    echo ${outvcf}
}

function run_gnomAD_genome_annotation {
    local region=$1
    local invcf=$2
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/b37-gnomAD-genome-exome-annotation-final.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting run_gnomAD_genome_annotation"
        echo ${outvcf}
        return 0;
    fi

    # annotation columns of interest
    local -a anno_cols=(
        'ID'
        'INFO/GNOMAD_GENOME_AC_AFR'
        'INFO/GNOMAD_GENOME_AC_AMR'
        'INFO/GNOMAD_GENOME_AC_ASJ'
        'INFO/GNOMAD_GENOME_AC_EAS'
        'INFO/GNOMAD_GENOME_AC_FIN'
        'INFO/GNOMAD_GENOME_AC_NFE'
        'INFO/GNOMAD_GENOME_AC_OTH'
        'INFO/GNOMAD_GENOME_AN_AFR'
        'INFO/GNOMAD_GENOME_AN_AMR'
        'INFO/GNOMAD_GENOME_AN_ASJ'
        'INFO/GNOMAD_GENOME_AN_EAS'
        'INFO/GNOMAD_GENOME_AN_FIN'
        'INFO/GNOMAD_GENOME_AN_NFE'
        'INFO/GNOMAD_GENOME_AN_OTH'
        'INFO/GNOMAD_GENOME_AC_raw'
        'INFO/GNOMAD_GENOME_AN_raw'
        'INFO/GNOMAD_GENOME_AC_POPMAX'
        'INFO/GNOMAD_GENOME_AN_POPMAX'
        'INFO/GNOMAD_GENOME_AS_RF'
        'INFO/GNOMAD_GENOME_AS_FilterStatus'
        'INFO/GNOMAD_GENOME_AS_RF_POSITIVE_TRAIN'
        'INFO/GNOMAD_GENOME_AS_RF_NEGATIVE_TRAIN'
    )
    local annotations=$(join_by , "${anno_cols[@]}")

    # the gnomAD genome annotation vcfs base location
    local base=/gscmnt/gc2802/halllab/gnomAD/release-170228/processed/post-vqsr-pipeline/genome

    # Figure out which gnomAD chromosomes vcfs are needed
    local -a gnomAD_chroms=($(${TABIX} --list-chroms ${invcf} | ${SORT} -N | ${UNIQ} | grep -v -P '^(GL|MT)'))
    log "[genome] gnomAD chromosomal vcfs to process: ${gnomAD_chroms[@]} ( items: ${#gnomAD_chroms[@]} )"

    local enter_vcf=${invcf}
    local exit_vcf=

    for chr in ${gnomAD_chroms[@]}; do
        log "Annotating with gnomAD genome chromosome: ${chr}"
        exit_vcf=${outdir}/b37-gnomAD-genome-annotation-c${chr}.vcf.gz

#        log "enter_vcf : ${enter_vcf}"
#        log "exit_vcf : ${exit_vcf}"

        # perform the chromosomal annotation

        if [[ "${chr}" == "Y" ]]; then
            # specially handle chromosome Y (there are no chr Y gnomAD annotations)
            log "gnomAD genome special chromosome Y handling"
            local cmd0="
            ${LN} -s ${enter_vcf} ${exit_vcf} && ${LN} -s ${enter_vcf}.tbi ${exit_vcf}.tbi
            "
            run_cmd "${cmd0}"
        else
            # determine the annotation vcf
            local anno_vcf=${base}/${chr}/
            anno_vcf+=c${chr}.gnomad.genomes.r2.0.1.sites.decompose.normalized.uniq.namespaced.vcf.gz

            if [[ ! -e ${anno_vcf} ]]; then
                die "Did not find annotation vcf: '${anno_vcf}'"
            fi

            local tmpvcf=${exit_vcf}.tmp
            local cmd1="
            ${BCFTOOLS} annotate \
                -a ${anno_vcf} \
                -c ${annotations} \
                -O z \
                -o ${tmpvcf} \
                ${enter_vcf} \
            "
            run_cmd "${cmd1}"
            tabix_and_finalize_vcf ${tmpvcf} ${exit_vcf}
        fi

        # setup for the next iteration
        enter_vcf=${exit_vcf}
    done

    log "final genome exit_vcf is: ${exit_vcf}"
    local cmd2="mv ${exit_vcf}.tbi ${outvcf}.tbi && mv ${exit_vcf} ${outvcf}"
    run_cmd "${cmd2}"

    echo ${outvcf}
}

function integrate_b37_annotations_to_b38 {
    local integrate_script=$1
    local b37_vcf=$2
    local b38_vcf=$3
    local annotation_type=$4

    local outdir=$(dirname ${b38_vcf})
    local outvcf=${outdir}/b38.genome.exome.nosamples.vcf.gz

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
        --update-id \
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
    local region=$1
    local b38_invcf=$2
    local b38_outvcf=$3
    local integrate_script=$4

    local scratch_dir=$(dirname ${b38_outvcf})/scratch
    mkdir -p ${scratch_dir}

    log "Remove samples on b38 input vcf"
    local b38_invcf_no_samples=$(prune_samples_on_b38_vcf ${b38_invcf} ${scratch_dir})
    log "Entering liftOver hg19"
    local hg19_vcf=$(run_liftover_hg19 ${b38_invcf_no_samples})
    log "Entering liftOver GRCh37"
    local grc37_vcf=$(run_liftover_grc37 ${hg19_vcf})
    log "Entering run_gnomAD_exome_annotation"
    local b37_gnomAD_exome_vcf=$(run_gnomAD_exome_annotation ${region} ${grc37_vcf})
    log "Entering run_gnomAD_genome_annotation"
    local b37_anno_vcf=$(run_gnomAD_genome_annotation ${region} ${b37_gnomAD_exome_vcf})
    log "Entering integrate b37 cadd annotations back to b38"
    local b38_anno_vcf=$(integrate_b37_annotations_to_b38 ${integrate_script} ${b37_anno_vcf} ${b38_invcf_no_samples} 'gnomAD')
    log "Add samples on b38 cadd annotated vcf"
    add_samples_on_b38_annotated_vcf ${b38_invcf} ${b38_anno_vcf} ${b38_outvcf}
}

function main {
    local chrom_region=$1
    local invcf=$2
    local outvcf=$3
    local integrate_script=$4

    if is_empty_vcf ${invcf} ; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${invcf} ${outvcf} ;
    else
        annotate_vcf ${chrom_region} ${invcf} ${outvcf} ${integrate_script};
    fi

    log 'All Done'
}

CHROM_REGION=$1
INVCF=$2
OUTVCF=$3
INTEGRATE_SCRIPT=$4

main ${CHROM_REGION} ${INVCF} ${OUTVCF} ${INTEGRATE_SCRIPT} ;
