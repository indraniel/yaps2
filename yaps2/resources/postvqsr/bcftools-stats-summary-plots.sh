#!/bin/bash

# Inputs:
#    INPUT_VCF_DIR     ($1) : the top level directory containing the per chromosome variant eval stats
#                             (e.g. $basedir/8-bcftools-stats)
#    OUTPUT_STATS_DIR  ($2) : the output directory that will contain the merged variant stats and plots
#                             (e.g. $basedir/8.1-bcftools-stats-summary)

# based on: March 2, 2016 3:40PM comment in JIRA issue BIO-1792
# https://jira.gsc.wustl.edu/browse/BIO-1792?focusedCommentId=108030&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-108030

# where to find a LaTeX distribution in a docker environment
if [ -e /.dockerenv ]; then
    echo "I am in a docker image"
    export PATH=/gscmnt/gc2802/halllab/idas/software/local/texlive-2015/2015/bin/x86_64-linux:$PATH
else
    echo "I am NOT in a docker image"
fi

# setup the python virtualenv environment (using python 2.7.10 with matplotlib==1.5.1)
BIO_1662=/gscmnt/gc2802/halllab/idas/jira/BIO-1662
source ${BIO_1662}/.env

# identify the input and output directories
INPUT_VCF_DIR=$1
OUTPUT_STATS_DIR=${2}/   # terminating '/' needed to create directory by plot-vcfstats

# the bcftools script that performs the work
PLOT_VCFSTATS=/gscuser/dlarson/src/bcftools/plot-vcfstats

# setup the desired perl
PERL=/gsc/scripts/opt/genome/current/user/bin/genome-perl

# get the relevant input stats files (generated by 7-bcftools-stats.sh)
INPUT_VCFS=${INPUT_VCF_DIR}/*/*.stats.out

# run the command
cmd="${PERL} ${PLOT_VCFSTATS} --prefix ${OUTPUT_STATS_DIR} ${INPUT_VCFS}"

echo $cmd
eval $cmd
