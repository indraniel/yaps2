#!/bin/bash

ZCAT=/bin/zcat
BCFTOOLS=/gsc/bin/bcftools1.2
BGZIP=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin/bgzip
TABIX=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin/tabix

BIO_1984=/gscmnt/gc2802/halllab/idas/jira/BIO-1984

# the 1000G annotation file
KGVCF=${BIO_1984}/data/manual/create-1000G-reformatted-af-annotations/
KGVCF+=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.decompose.normalize.reheader.w_ids.reformatted_pop_af.vcf.gz

INVCF=$1
OUTVCF=$2

function is_empty_vcf {
    local count=$(${ZCAT} ${INVCF} | head -n 1000 | wc -l)
    if [[ "${count}" -gt "0" ]]; then
        return 1
    else
        return 0
    fi
}

function copy_over_vcf_headers {
    local region=$(basename $(dirname ${INVCF}))
    local project_path=$(dirname $(dirname $(dirname ${INVCF})))
    local base_vcf=${project_path}/3-decompose-normalize-uniq/${region}/combined.c${region}.vcf.gz
    ${BCFTOOLS} view -h ${base_vcf} | ${BGZIP} -c >${OUTVCF}
    ${TABIX} ${OUTVCF}
}

function annotate_vcf {
    ${BCFTOOLS} annotate \
        -a ${KGVCF} \
        -c ID,INFO/1KG_EAS_AF,INFO/1KG_EUR_AF,INFO/1KG_AFR_AF,INFO/1KG_AMR_AF,INFO/1KG_SAS_AF \
        -O z \
        -o ${OUTVCF} \
        ${INVCF} \
        && ${TABIX} -p vcf -f ${OUTVCF}
}

function main {
    if is_empty_vcf ; then
        echo "No variants to process copying headers over..."
        copy_over_vcf_headers ;
    else
        annotate_vcf ;
    fi
}

main ;
