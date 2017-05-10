#!/bin/bash

set -eo pipefail

AWK=/usr/bin/awk

JAVA=/gapp/x64linux/opt/java/jdk/jdk1.8.0_60/bin/java
PICARD=/gscmnt/gc2802/halllab/idas/software/picard/picard.2.9.0.jar

PYTHON=$(which python) # if run inside yaps2 pipeline, then should be getting the virtualenv python

BGZIP=/gscmnt/gc2802/halllab/idas/software/local/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/local/bin/tabix
BCFTOOLS=/gscmnt/gc2802/halllab/idas/software/local/bcftools1.4

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

# function subtract_vcf_samples {
#     local invcf=$1
#     local new_vcf_name=$2
#     local outdir=$(dirname ${invcf})
#     local outvcf=${outdir}/${new_vcf_name}
#     local tmpvcf=$outvcf.tmp
# 
#     if [[ -e "${outvcf}" ]]; then
#         log "shortcutting subtract_vcf_samples (${new_vcf_name})"
#         echo ${outvcf}
#         return 0;
#     fi
# 
#     # add the header to the new test vcf (except for the last header line)
#     cmd1="${BCFTOOLS} view -h ${invcf} | head -n -1 > ${tmpvcf}"
#     run_cmd "${cmd1}"
# 
#     # correctly add the last header line
#     cmd2="${BCFTOOLS} view -h ${invcf} | tail -n 1 | cut -f 1-8 >>${tmpvcf}"
#     run_cmd "${cmd2}"
# 
#     # add the first few samples / variants
#     cmd3="${BCFTOOLS} view -H ${invcf} | cut -f 1-8 >>${tmpvcf}"
#     run_cmd "${cmd3}"
# 
#     # tabix & bgzip the file
#     local cmd4="
#     ${TABIX} -p vcf -f ${tmpvcf} \
#         && mv ${tmpvcf}.tbi ${outvcf}.tbi \
#         && mv ${tmpvcf} ${outvcf}
#     "
#     run_cmd "${cmd4}"
#     echo ${outvcf}
# }

function vcf_subtract_samples {
    local vcf=$1

    cat <(${BCFTOOLS} view -h ${vcf} | head -n -1) \
        <(${BCFTOOLS} view -h ${vcf} | tail -n 1 | cut -f1-8) \
        <(${BCFTOOLS} view -H ${vcf} | cut -f1-8)
}

function vcf_add_samples {
    local no_samples_vcf=$1
    local samples_vcf=$2

    cat <(${BCFTOOLS} view -h ${no_samples_vcf} | head -n -1) \
        <(${BCFTOOLS} view -h ${samples_vcf} | tail -n 1 ) \
        <(paste <(${BCFTOOLS} view -H ${no_samples_vcf}) \
                <(${BCFTOOLS} view -H ${samples_vcf} | cut -f9-))
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

function prune_grc37vcf_samples {
    local invcf=$1
    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/grc37.nosample.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting prune_grc37vcf_samples"
        echo ${outvcf}
        return 0;
    fi

    local tmptsv=${outvcf}.tmp

	local cmd1="cat <(vcf_subtract_samples ${invcf}) | ${BGZIP} -c > ${tmpvcf}"
	run_cmd "${cmd1}"

	tabix_and_finalize_vcf ${tmpvcf} ${outvcf}

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

    local perl_path=/gscmnt/gc2802/halllab/idas/software/std-perl
    local htslib_path=/gscmnt/gc2802/halllab/idas/software/local/bin
    local awk_path=/gscmnt/gc2802/halllab/idas/software/std-awk
    export PERL5LIB=/gscmnt/gc2719/halllab/src/speedseq/bin:${PERL5LIB}
    export PATH=${perl_path}:${awk_path}:${htslib_path}:${PATH}

    local cadd_dir=/gscmnt/gc2802/halllab/abelhj/CADD/CADD_v1.2
    local cadd_score_script=${cadd_dir}/bin/score.sh

    local cmd="/bin/bash ${cadd_score_script} ${invcf} ${tmptsv} && mv ${tmptsv} ${tsv}"
    run_cmd "${cmd}"

    echo ${tsv}
}

function paste_cadd {
    local merge_cadd_script=$1
    local invcf=$2
    local caddtsv=$3
    local fullvcf=$4

    local outdir=$(dirname ${invcf})
    local outvcf=${outdir}/grc37.vep.cadd.vcf.gz

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting paste_cadd"
        echo ${outvcf}
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # core of "paste_cadd.sh"
# old approach    
#    local cmd1="
#    cat <(${PYTHON} ${merge_cadd_script} -H <(zcat ${invcf}) <(zcat ${caddtsv})) \
#        <(paste <(${PYTHON} ${merge_cadd_script} --no-header <(zcat ${invcf}) <(zcat ${caddtsv})) \
#                <(zcat ${fullvcf} | ${AWK} '{if(\$5!=\".\") print \$0}' | grep -v '^##' | cut -f9- )) \
#        | ${BGZIP} -c \
#        > ${tmpvcf}
#    "

# # new approach (based on changes in gatk-workflow.git:scripts/paste_cadd.sh)
#     local cmd1="
#     cat <(python2.7 ${merge_cadd_script} -H <(zcat ${invcf}) <(zcat ${caddtsv})) \
#         <(paste <(python2.7 ${merge_cadd_script} --no-header \
#                     <(cat <(zcat ${invcf} | grep \"^#\") \
#                           <(zcat ${invcf} \
#                               | grep -v \"^#\" \
#                               | sort -k1,1 -k2,2n -k4,4 -k5,5)) \
#                     <(zcat ${caddtsv})) \
#                 <(zcat ${fullvcf} \
#                     | grep -v \"^##\" \
#                     | sort -k1,1 -k2,2n -k4,4 -k5,5 \
#                     | ${AWK} '{if(\$5!=\".\") print \$0}' \
#                     | cut -f9- )) \
#     | ${BGZIP} -c \
#     > ${tmpvcf}
#     "
#     run_cmd "${cmd1}"

# even newer approach (defer merging of the samples for a later step)
    local cmd1="
    cat <(python2.7 ${merge_cadd_script} -H <(zcat ${invcf}) <(zcat ${caddtsv})) \
        <(python2.7 ${merge_cadd_script} --no-header \
                    <(cat <(zcat ${invcf} | grep \"^#\") \
                          <(zcat ${invcf} \
                              | grep -v \"^#\" \
                              | sort -k1,1 -k2,2n -k4,4 -k5,5)) \
                    <(zcat ${caddtsv})) \
    | ${BGZIP} -c \
    > ${tmpvcf}
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

    ${PYTHON} ${integrate_script} \
        --b38-vcf=$(vcf_subtract_samples ${b38_vcf})  \
        --annotated-b37-vcf=${b37_vcf} \
        --auto-fill \
        --annotation-type=${annotation_type}
}

function migrate_b37_annotations_to_b38 {
	local original_b38_input_vcf=$1
    local final_b38_annotated_output_vcf=$2
    local b37vcf=$3
    local integrate_script=$4

    if [[ -e "${final_b38_annotated_output_vcf}" ]]; then
        log "shortcutting migrate_b37_annotations_to_b38"
        return 0;
    fi

    local tmpvcf=${final_b38_annotated_output_vcf}.tmp

    local cmd1="
	vcf_add_samples \
		<(integrate_b37_annotations_to_b38 \
             ${integrate_script} \
             ${b37vcf} \
		     ${original_b38_input_vcf} \
             'cadd') \
        ${original_b38_input_vcf} \
    | ${BGZIP} -c \
    > ${tmpvcf}
    "
    run_cmd "${cmd1}"

	tabix_and_finalize_vcf ${tmpvcf} ${final_b38_annotated_output_vcf}
}

function annotate_vcf {
    local invcf=$1
    local outvcf=$2
    local merge_script=$3
    local migrate_script=$4

    log "Entering liftOver hg19"
    local hg19vcf=$(run_liftover_hg19 ${invcf})
    log "Entering liftOver GRCh37"
    local grc37vcf=$(run_liftover_grc37 ${hg19vcf})
    log "Remove samples on GRCh37 liftOver vcf"
    local no_sample_grc37_vcf=$(prune_grc37vcf_samples ${grc37vcf})
    log "Entering run_cadd"
    local tsv=$(run_cadd ${no_sample_grc37vcf})
    log "Entering paste_cadd"
    local b37vcf=$(paste_cadd ${merge_script} ${no_sample_grc37_vcf} ${tsv} ${grc37vcf} )
    log "Entering migrate b37 annotations back to b38"
    migrate_b37_annotations_to_b38 ${invcf} ${outvcf} ${b37vcf} ${migrate_script}
}

function main {
    local invcf=$1
    local outvcf=$2
    local merge_script=$3
    local migrate_script=$4

    if is_empty_vcf ${invcf} ; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${invcf} ${outvcf} ;
    else
        annotate_vcf ${invcf} ${outvcf} ${merge_script} ${migrate_script};
    fi

    log 'All Done'
}

INVCF=$1
OUTVCF=$2
MERGE_SCRIPT=$3
MIGRATE_B37_ANNOTATIONS_TO_B38_SCRIPT=$4

main ${INVCF} ${OUTVCF} ${MERGE_SCRIPT} ${MIGRATE_B37_ANNOTATIONS_TO_B38_SCRIPT};
