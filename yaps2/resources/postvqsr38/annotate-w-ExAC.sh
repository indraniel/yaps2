#!/bin/bash

BCFTOOLS=/gsc/bin/bcftools1.2
TABIX=/gscmnt/gc2802/halllab/idas/software/vep/local/htslib-1.3.2/bin/tabix

BIO_1984=/gscmnt/gc2802/halllab/idas/jira/BIO-1984

# the ExAC annotation file
KGVCF=${BIO_1984}/data/manual/create-ExAC-reformatted-af-annotations/
KGVCF+=ExAC.r0.3.1.sites.simpleAF.multiallelic.reformatted_pop_af.vcf.gz

INVCF=$1
OUTVCF=$2

function is_empty_vcf {
    local count=$(${BCFTOOLS} view -H ${INVCF} | head -n 1000 | wc -l)
    if [[ "${count}" -gt "0" ]]; then
        return 1
    else
        return 0
    fi
}

function copy_over_vcf {
    cp -v ${INVCF} ${OUTVCF}
    cp -v ${INVCF}.tbi ${OUTVCF}.tbi
}

function annotate_vcf {
    ${BCFTOOLS} annotate \
        -a ${KGVCF} \
        -c INFO/EXAC_AFR_AF,INFO/EXAC_AMR_AF,INFO/EXAC_Adj_AF,INFO/EXAC_EAS_AF,INFO/EXAC_FIN_AF,INFO/EXAC_NFE_AF,INFO/EXAC_OTH_AF,INFO/EXAC_SAS_AF \
        -O z \
        -o ${OUTVCF} \
        ${INVCF} \
        && ${TABIX} -p vcf -f ${OUTVCF}
}

function main {
    if is_empty_vcf ; then
        echo "No variants to process copying files over..."
        copy_over_vcf ;
    else
        annotate_vcf ;
    fi
}

main ;
