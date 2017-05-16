.PHONY: clean

SHELL := /bin/bash
MKFILE_PATH := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# assume we're running inside an LSF job
WORKSPACE := /tmp/$(LSB_JOBID).tmpdir/1kg-test

# location of a python 2.7.x virtualenv 
export VIRTUAL_ENV := /gscmnt/gc2802/halllab/idas/laboratory/yaps2-cadd-vep-test/test-venv

PRG_DIR          := $(MKFILE_PATH)/../../yaps2/resources/postvqsr38
SCRIPT           := $(PRG_DIR)/annotate-allele-balances.py
ORIG_SCRIPT      := /gscmnt/gc2802/halllab/abelhj/cvd0223/qc/add_ab_both.0515.pl
BCFTOOLS         := /gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4

ORIG_INPUT_VCF   := /gscmnt/gc2802/halllab/aregier/jira/BIO-2228/decomposed.vcf.gz
TEST_INPUT_VCF   := $(WORKSPACE)/in/decomposed.vcf.gz
TEST_OUTPUT_VCF  := $(WORKSPACE)/out/b38.ab.annotated.vcf
TEST_ORIG_OUTPUT_VCF  := $(WORKSPACE)/out/b38.ab.annotated.orig.vcf

TEST_OUTPUT_STATS := $(WORKSPACE)/out/b38.ab.stats
TEST_ORIG_STATS := $(WORKSPACE)/out/b38.ab.orig.stats
COMPARE_FILE := $(WORKSPACE)/out/compare.dat

run:
	python $(SCRIPT) $(TEST_INPUT_VCF) >$(TEST_OUTPUT_VCF)
	zcat $(TEST_INPUT_VCF) | perl $(ORIG_SCRIPT) >$(TEST_ORIG_OUTPUT_VCF)
	$(BCFTOOLS) query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%INFO/HetAB\t%INFO/HetHomAltAB\n' $(TEST_OUTPUT_VCF) > $(TEST_OUTPUT_STATS)
	$(BCFTOOLS) query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%INFO/AB\t%INFO/AB_HOM\n' $(TEST_ORIG_OUTPUT_VCF) > $(TEST_ORIG_STATS)
	paste $(TEST_ORIG_STATS) $(TEST_OUTPUT_STATS) >$(COMPARE_FILE)

setup:
	mkdir -p $(WORKSPACE)/{in,out}
	cp $(ORIG_INPUT_VCF) $(WORKSPACE)/in

clean:
	rm -rf $(WORKSPACE)
