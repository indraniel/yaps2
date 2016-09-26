.PHONY: clean reformat check-env create-project-dir

# set the shell to use bash
SHELL := /bin/bash

# needed inputs
INPUT_VCF :=                # aka "trios.vcf.gz"
PRJ_DIR  := 
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
empty_vcf := $(PRJ_DIR)/empty-vcf

INPUT_VCF_VARIANT_COUNT := $(shell $(BCFTOOLS) view $(INPUT_VCF) | grep -v '^\#' | head -n 1 | wc -l )

ifeq ($(shell test $(INPUT_VCF_VARIANT_COUNT) -ge 1; echo $$?),0)
EMPTY_INPUT_VCF_FLAG :=
else
EMPTY_INPUT_VCF_FLAG := "empty"
endif

all: check-env create-project-dir $(unfiltered-mie-var-txt)

$(unfiltered-mie-var-txt): $(unfiltered-info-variants-var) $(unfiltered-fmendel)
ifndef EMPTY_INPUT_VCF_FLAG
	$(BASH) $(MAKE_TABLE) $(unfiltered-info-variants-var) $(unfiltered-fmendel) > $(unfiltered-mie-var-txt)
else
	echo "$(INPUT_VCF) has no variants!"
	touch $(empty_vcf)
endif

$(unfiltered-info-variants-var): $(TRIO_FAM) $(INPUT_VCF)
	$(BASH) $(SUMMARIZE_INFORMATIVE_VAR_ONLY) $(TRIO_FAM) ${INPUT_VCF} > $(unfiltered-info-variants-var)

$(unfiltered-mendel) $(unfiltered-fmendel) $(unfiltered-imendel) $(unfiltered-lmendel): $(unfiltered-bed) $(unfiltered-bim) $(unfiltered-fam)
ifndef EMPTY_INPUT_VCF_FLAG
	# excluding anything that isn't on an autosome as I suspect XY etc are not well handled
	$(PLINK) --bfile $(unfiltered) --mendel --allow-extra-chr --out $(unfiltered) --chr 1-22
	$(eval PLINKOUT := $(unfiltered-mendel) $(unfiltered-fmendel) $(unfiltered-imendel) $(unfiltered-lmendel))
	$(call reformat,$(PLINKOUT))
else
	touch $(unfiltered-fmendel)
endif

$(unfiltered-bed) $(unfiltered-bim) $(unfiltered-fam): $(INPUT_VCF) $(TRIO_FAM)
ifndef EMPTY_INPUT_VCF_FLAG
	# allow-extra-chr is there to allow unplaced contigs
	# double id is a reasonable default, but we will replace the ouput .fam file with our own later below
	$(PLINK) --vcf $(INPUT_VCF) --double-id --make-bed --out $(unfiltered) --allow-extra-chr
	$(RM) $(unfiltered-nosex)
	$(CP) $(TRIO_FAM) $(unfiltered-fam)
else
	touch $(unfiltered-bed) $(unfiltered-bim) $(unfiltered-fam)
	$(CP) $(TRIO_FAM) $(unfiltered-fam)
endif

define reformat
	# This takes each plink output file and translates it to use tabs.
	# Commands taken from https://www.cog-genomics.org/plink2/other#tabspace
	$(foreach file,$1,cat $(file) \
		| sed 's/^[[:space:]]*//g' \
		| sed 's/[[:space:]]*$$//g' \
		| tr -s ' ' '\t' > $(file).tmp && mv $(file).tmp $(file);)
endef

create-project-dir:
	mkdir -p $(PRJ_DIR)

clean: check-env
	$(RM) $(PRJ_DIR)

check-env:
ifndef PRJ_DIR
	$(error PRJ_DIR is undefined!)
endif

ifndef TRIO_FAM
	$(error TRIO_FAM is undefined!)
endif

ifndef INPUT_VCF
	$(error INPUT_VCF is undefined!)
endif
