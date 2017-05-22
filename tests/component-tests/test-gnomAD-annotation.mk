.PHONY: clean

SHELL := /bin/bash
MKFILE_PATH := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# assume we're running inside an LSF job
WORKSPACE := /tmp/$(LSB_JOBID).tmpdir/gnomAD-test

# location of a python 2.7.x virtualenv 
export VIRTUAL_ENV := /gscmnt/gc2802/halllab/idas/laboratory/yaps2-cadd-vep-test/test-venv

PRG_DIR           := $(MKFILE_PATH)/../../yaps2/resources/postvqsr38
SCRIPT            := $(PRG_DIR)/annotate-w-gnomAD.sh
REGION            := "chr10:41007921-49725146"
INTEGRATE_SCRIPT  := $(PRG_DIR)/integrate-b37-annotations-to-b38.py

ORIG_INPUT_VCF   := /gscmnt/gc2802/halllab/aregier/jira/BIO-2228/decomposed.vcf.gz
TEST_INPUT_VCF   := $(WORKSPACE)/in/decomposed.vcf.gz
TEST_OUTPUT_VCF  := $(WORKSPACE)/out/b38.final.gnomAD.annotated.vcf.gz

run:
	source $(VIRTUAL_ENV)/bin/activate \
		&& bash $(SCRIPT) $(REGION) $(TEST_INPUT_VCF) $(TEST_OUTPUT_VCF) $(INTEGRATE_SCRIPT)

setup:
	mkdir -p $(WORKSPACE)/{in,out}
	cp $(ORIG_INPUT_VCF) $(WORKSPACE)/in

clean:
	rm -rf $(WORKSPACE)
