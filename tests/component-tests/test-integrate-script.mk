.PHONY: clean

SHELL := /bin/bash
PYTHON := /opt/pyenv/.pyenv/versions/2.7.10/bin/python
MKFILE_PATH := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# assume we're running inside an LSF job
WORKSPACE := /tmp/$(LSB_JOBID).tmpdir/integrate-id-test

# location of a python 2.7.x virtualenv 
export VIRTUAL_ENV := /gscmnt/gc2802/halllab/idas/laboratory/yaps2-cadd-vep-test/test-venv

PRG_DIR           := $(MKFILE_PATH)/../../yaps2/resources/postvqsr38
INTEGRATE_SCRIPT  := $(PRG_DIR)/integrate-b37-annotations-to-b38.py

ORIG_B37_INPUT_VCF   := $(MKFILE_PATH)/data/integrate.id.test.grc37.cadd.vcf.gz
ORIG_B38_INPUT_VCF   := $(MKFILE_PATH)/data/integrate.id.test.b38.nosamples.vcf.gz
TEST_B37_INPUT_VCF   := $(WORKSPACE)/in/integrate.id.test.grc37.cadd.vcf.gz
TEST_B38_INPUT_VCF   := $(WORKSPACE)/in/integrate.id.test.b38.nosamples.vcf.gz
TEST_OUTPUT_VCF      := $(WORKSPACE)/out/b38.final.annotated.vcf.gz

run:
	$(PYTHON) $(INTEGRATE_SCRIPT) \
		--b38-vcf=$(TEST_B38_INPUT_VCF) \
		--annotated-b37-vcf=$(TEST_B37_INPUT_VCF) \
		--auto-fill \
		--annotation-type=cadd \
		--update-id \
		> $(TEST_OUTPUT_VCF)

setup:
	mkdir -p $(WORKSPACE)/{in,out}
	cp $(ORIG_B37_INPUT_VCF) $(WORKSPACE)/in
	cp $(ORIG_B37_INPUT_VCF).tbi $(WORKSPACE)/in
	cp $(ORIG_B38_INPUT_VCF) $(WORKSPACE)/in
	cp $(ORIG_B38_INPUT_VCF).tbi $(WORKSPACE)/in

clean:
	rm -rf $(WORKSPACE)
