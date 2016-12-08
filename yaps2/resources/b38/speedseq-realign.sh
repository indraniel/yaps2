#!/gsc/bin/bash

# based on: /gscmnt/gc2802/halllab/aregier/jira/BIO-1875/scripts/run_pipeline_bams.sh

SPEEDSEQ=/gscmnt/sata849/info/speedseq_testing/v0.2.0-gms-testing/speedseq/bin/speedseq

create_dir() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        echo "Creating directory: $dir"
        mkdir -p "$dir"
    fi
}

setup_virtualenv() {
    if [ $VIRTUAL_ENV ]; then
        echo "Found ${VIRTUAL_ENV} -- activating"
        source ${VIRTUAL_ENV}/bin/activate
    fi
}

OUTPUT_PREFIX=$1
TEMP=$2
REF=$3
BAM_STRING="${@:4}"

echo "BAM_STRING:$BAM_STRING"
echo "REF:$REF"

create_dir $TEMP
setup_virtualenv

cmd="${SPEEDSEQ} realign -o ${OUTPUT_PREFIX} -t 8 -T ${TEMP} ${REF} ${BAM_STRING}"
echo "EXECUTING: ${cmd}"
eval ${cmd}
