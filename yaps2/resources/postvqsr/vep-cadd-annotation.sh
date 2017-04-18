#!/bin/bash

# based on: /gscmnt/gc2802/halllab/aregier/gatk-workflow/scripts/run_cadd_wrapper.sh
#           /gscmnt/gc2802/halllab/aregier/gatk-workflow/scripts/run_score.sh
#           /gscmnt/gc2802/halllab/aregier/gatk-workflow/scripts/paste_cadd_wrapper.sh
#           /gscmnt/gc2802/halllab/aregier/gatk-workflow/scripts/paste_cadd.sh

# http://stackoverflow.com/questions/9893667/is-there-a-way-to-write-a-bash-function-which-aborts-the-whole-execution-no-mat
trap "exit 1" TERM
export TOP_PID=$$

AWK=/usr/bin/awk
PYTHON=$(which python) # if run inside yaps2 pipeline, then should be getting the virtualenv python
BGZIP=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin/tabix
BCFTOOLS=/gscmnt/gc2719/halllab/bin/bcftools

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
        kill -s TERM $TOP_PID
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

    local vep=/gscmnt/gc2719/halllab/src/speedseq/bin/variant_effect_predictor.pl
    local vep_cache=/gscmnt/gc2719/halllab/src/speedseq/annotations/vep_cache
    local vep_plugin_dir=/gscmnt/gc2802/aregier/bin/loftee
    local vep_plugin_args="LoF,human_ancestor_fa:/gscmnt/gc2802/aregier/bin/loftee/human_ancestor.fa.gz,conservation_file:/gscmnt/gc2802/aregier/bin/loftee/phylocsf.sql"
    local vep_fields="Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,LoF,LoF_filter,LoF_flags,LoF_info"

    local perl_path=/gscmnt/gc2802/halllab/idas/software/std-perl
    local htslib_path=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin
    export PERL5LIB=/gscmnt/gc2719/halllab/src/speedseq/bin:${PERL5LIB}
    export PATH=${perl_path}:${htslib_path}:${PATH}

    local cmd1="
    zcat ${invcf} \
        | cut -f 1-8 \
        | ${vep} \
            --force_overwrite \
            --offline \
            --fork 12 \
            --cache \
            --dir_cache ${vep_cache} \
            --dir_plugins ${vep_plugin_dir} \
            --plugin ${vep_plugin_args} \
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

function run_cadd {
    local invcf=$1
    local outdir=$2

    local tsv=${outdir}/cadd-annotation.tsv.gz

    if [[ -e "${tsv}" ]]; then
        log "shortcutting run_cadd"
        echo ${tsv}
        return 0;
    fi

    local tmptsv=${tsv}.tmp

    local perl_path=/gscmnt/gc2802/halllab/idas/software/std-perl
    local htslib_path=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin
    local awk_path=/gscmnt/gc2802/halllab/idas/software/std-awk
    export PERL5LIB=/gscmnt/gc2719/halllab/src/speedseq/bin:${PERL5LIB}
    export PATH=${perl_path}:${awk_path}:${htslib_path}:${PATH}

    local cadd_dir=/gscmnt/gc2802/halllab/abelhj/CADD/CADD_v1.2
    cadd_score_script=${cadd_dir}/bin/score.sh

    local cmd="/bin/bash ${cadd_score_script} ${invcf} ${tmptsv} && mv ${tmptsv} ${tsv}"
    run_cmd "${cmd}"

    echo ${tsv}
}

function paste_cadd {
    local tsvfile=$1
    local vepvcf=$2
    local origvcf=$3
    local outvcf=$4
    local merge_cadd_script=$5

    if [[ -e "${outvcf}" ]]; then
        log "shortcutting paste_cadd"
        return 0;
    fi

    local tmpvcf=${outvcf}.tmp

    # core of "paste_cadd.sh"
    local cmd1="
    cat <(${PYTHON} ${merge_cadd_script} -H <(zcat ${vepvcf}) <(zcat ${tsvfile})) \
        <(paste <(${PYTHON} ${merge_cadd_script} --no-header <(zcat ${vepvcf}) <(zcat ${tsvfile})) \
                <(zcat ${origvcf} | ${AWK} '{if(\$5!=\".\") print \$0}' | cut -f9- | grep -v '^#')) \
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
}

function annotate_vcf {
    local invcf=$1
    local outvcf=$2
    local merge_script=$3

    local outdir=$(dirname ${outvcf})

    log "Entering run_vep"
    local vepvcf=$(run_vep ${invcf} ${outdir})
    log "Entering run_cadd"
    local tsv=$(run_cadd ${vepvcf} ${outdir})
    log "Entering paste_cadd"
    paste_cadd ${tsv} ${vepvcf} ${invcf} ${outvcf} ${merge_script}
}

function main {
    local invcf=$1
    local outvcf=$2
    local merge_script=$3

    if is_empty_vcf ${invcf} ; then
        log "No variants to process. Copying files over..."
        copy_over_vcf ${invcf} ${outvcf} ;
    else
        annotate_vcf ${invcf} ${outvcf} ${merge_script};
    fi

    log 'All Done'
}

INVCF=$1
OUTVCF=$2
MERGE_SCRIPT=$3

main ${INVCF} ${OUTVCF} ${MERGE_SCRIPT};
