.PHONY: clean reformat check-env

# needed inputs
INPUT_VCF :=                # aka "trios.vcf.gz"
PRJ_DIR  := 
COHORT_SEX :=               # aka /gscmnt/gc2802/halllab/dlarson/svtools_tests/benchmarking_with_trios/trios_plus_finns/data/reclass/cohort.sex
# Note that the order of this .fam (trio.fam) file HAS TO MATCH the order of the samples in the VCF
TRIO_FAM :=                 # aka "cohort.fam"

# relevant directories
PRG_ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

# programs
PLINK                          := /gscuser/dlarson/bin/plink
BCFTOOLS                       := /gsc/bin/bcftools1.2
MAKE_TABLE                     := $(PRG_ROOT_DIR)/make_table.sh
BASH                           := /bin/bash
RM                             := /bin/rm -rf
CP                             := /bin/cp -f
SUMMARIZE_INFORMATIVE_VAR_ONLY := $(PRG_ROOT_DIR)/summarize_informative_var_only.sh

# files of interest
unfiltered := $(PRJ_DIR)/unfiltered
unfiltered-nosex := $(unfiltered).nosex
unfiltered-fam := $(unfiltered).fam
unfiltered-bed := $(unfiltered).bed
unfiltered-bim := $(unfiltered).bim
unfiltered-mie-var-txt := $(unfiltered).mie.var.txt
unfiltered-info-variants-var := $(unfiltered).info_variants.var
unfiltered-mendel := $(unfiltered).mendel
unfiltered-fmendel := $(unfiltered).fmendel
unfiltered-imendel := $(unfiltered).imendel
unfiltered-lmendel := $(unfiltered).lmendel

all: check-env $(unfiltered-mie-var-txt) $(unfiltered-info-variants-var)

$(unfiltered-mie-var-txt): $(unfiltered-info-variants-var) $(unfiltered-fmendel)
	$(BASH) $(MAKE_TABLE) $(unfiltered-info-variants-var) $(unfiltered-fmendel) > $(unfiltered-mie-var-txt)

$(unfiltered-info-variants-var): $(TRIO_FAM) $(INPUT_VCF)
	$(BASH) $(SUMMARIZE_INFORMATIVE_VAR_ONLY) $(TRIO_FAM) ${INPUT_VCF} > $(unfiltered-info-variants-var)

$(unfiltered-mendel) $(unfiltered-fmendel) $(unfiltered-imendel) $(unfiltered-lmendel): $(unfiltered-bed) $(unfiltered-bim) $(unfiltered-fam)
	# excluding anything that isn't on an autosome as I suspect XY etc are not well handled
	$(PLINK) --bfile $(unfiltered) --mendel --allow-extra-chr --out $(unfiltered) --chr 1-22

$(unfiltered-bed) $(unfiltered-bim) $(unfiltered-fam): $(INPUT_VCF) $(TRIO_FAM)
	# allow-extra-chr is there to allow unplaced contigs
	# double id is a reasonable default, but we will replace the ouput .fam file with our own later below
	$(PLINK) --vcf $(INPUT_VCF) --double-id --make-bed --out $(unfiltered) --allow-extra-chr
	$(RM) $(unfiltered-nosex)
	$(CP) $(TRIO_FAM) $(unfiltered-fam)

PLINKOUT = $(unfiltered-mendel) $(unfiltered-fmendel) $(unfiltered-imendel) $(unfiltered-lmendel)
reformat:
	# This takes each plink output file and translates it to use tabs.
	# Commands taken from https://www.cog-genomics.org/plink2/other#tabspace
	$(foreach file,$(PLINKOUT),cat $(file) \
		| sed 's/^[[:space:]]*//g' \
		| sed 's/[[:space:]]*$$//g' \
		| tr -s ' ' '\t' > $(file).tmp && mv $(file).tmp $(file);)

clean:
	$(RM) $(PRJ_DIR)

check-env:
ifndef PRJ_DIR
	$(error PRJ_DIR is undefined!)
endif

ifndef COHORT_SEX
	$(error COHORT_SEX is undefined!)
endif

ifndef TRIO_FAM
	$(error TRIO_FAM is undefined!)
endif

ifndef INPUT_VCF
	$(error INPUT_VCF is undefined!)
endif
