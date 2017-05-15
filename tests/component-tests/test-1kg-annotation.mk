.PHONY: clean

SHELL := /bin/bash
MKFILE_PATH := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# assume we're running inside and LSF job
WORKSPACE := /tmp/$(LSB_JOBID).tmpdir/1kg-test

# location of a python 2.7.x virtualenv 
VIRTUAL_ENV := /gscmnt/gc2802/halllab/idas/laboratory/yaps2-cadd-vep-test/test-venv

PRG_DIR          := $(MKFILE_PATH)/../../yaps2/resources/postvqsr38
SCRIPT           := $(PRG_DIR)/annotate-w-1000G.sh
INTEGRATE_SCRIPT := $(PRG_DIR)/integrate-b37-annotations-to-b38.py

ORIG_INPUT_VCF   := /gscmnt/gc2802/halllab/aregier/jira/BIO-2228/decomposed.vcf.gz
TEST_INPUT_VCF   := $(WORKSPACE)/in/decomposed.vcf.gz
TEST_OUTPUT_VCF  := $(WORKSPACE)/out/b38.final.1k.annotated.vcf.gz

run:
	bash $(SCRIPT) $(TEST_INPUT_VCF) $(TEST_OUTPUT_VCF) $(INTEGRATE_SCRIPT)

setup:
	mkdir -p $(WORKSPACE)/{in,out}
	cp $(ORIG_INPUT_VCF) $(WORKSPACE)/in

clean:
	rm -rf $(WORKSPACE)
