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

function die {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo -e "[ ${timestamp} ] ERROR: $@" >&2
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

function usage {
    msg="Usage: $0 -m /path/to/final-merged-file.vcf.gz <one or more /path/to/vcf-files-to-merge.vcf.gz> ..."
    echo "${msg}";
}

function parse_args {
    local final_merged_vcf;
    while getopts "m:" arg; do
        case "${arg}" in
            m)
                final_merged_vcf=${OPTARG}
                ;;
            *)
                local msg=$(usage)
                die ${msg}
                ;;
        esac
    done
    shift $((OPTIND - 1))

    if [ -z ${final_merged_vcf} ]; then
        local msg=$(usage)
        die "[err] Did not find the argument to the -m parameter:\n${msg}"
    fi
    log "final_merged_vcf is: '${final_merged_vcf}'"
    log "input vcfs are: '${@}'"
    local count=1
    for f in ${@}; do
        log "${count} : ${f}"
        count=$((count + 1))
    done

    # returning bash associative array from function trick:
    # http://notes-matthewlmcclure.blogspot.com/2009/12/return-array-from-bash-function-v-2.html
    local -A params=( [final_merged_vcf]=${final_merged_vcf} [input_vcfs]="${@}" )
    declare -p params | sed -e 's/^declare -A params=//'
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

function merge {
    local params_string=$@
    eval "local -A params=${params_string}"
    local out_vcf=${params[final_merged_vcf]}
    local tmp_vcf=${out_vcf}.tmp
    local input_vcfs=${params[input_vcfs]}

    if [[ -e "${out_vcf}" ]]; then
        log "merged vcf already exists -- shortcutting merging"
        return 0;
    fi

    if [[ ${#input_vcfs} == 0 ]]; then
        die "[err] Did not find any input vcfs to merge!"
        exit 1
    fi

    log "(merge) calling bcftools concat"
    local cmd="
    ${BCFTOOLS} concat \
        -a ${input_vcfs} \
        -O z \
        -o ${tmp_vcf}
    "
    run_cmd "${cmd}"
    log "(merge) tabix and finalizing vcf"
	tabix_and_finalize_vcf ${tmp_vcf} ${out_vcf}

    echo ${out_vcf}
}

function main {
    local params_string=$(parse_args "$@")
    echo "params_string is: ${params_string}"
    local merged_vcf=$(merge ${params_string})
    log "Got the merged vcf: ${merged_vcf}"
}

main "$@"
