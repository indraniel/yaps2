.PHONY: clean check-env

# set the shell to use bash
SHELL := /bin/bash

# update the environment to include directory of custom install programs
EXTERNAL_SOFTWARE := /gscmnt/gc2802/halllab/idas/software/local
export PATH := $(EXTERNAL_SOFTWARE)/bin:$(PATH)
export LD_LIBRARY_PATH := $(EXTERNAL_SOFTWARE)/lib:$(LD_LIBRARY_PATH)

# needed inputs
PRJ_DIR  := 
INPUT_PED :=
INPUT_MAP :=

# programs
PERL := /usr/bin/perl
CONVERTF := $(EXTERNAL_SOFTWARE)/bin/convertf
SMARTPCA := $(EXTERNAL_SOFTWARE)/bin/smartpca.perl

# files of interest
merge-par := $(PRJ_DIR)/merged.par
genotype-output := $(PRJ_DIR)/merged.eigenstrat.geno
snp-output := $(PRJ_DIR)/merged.eigenstrat.snp
indiv-output := $(PRJ_DIR)/merged.eigenstrat.indiv
smartpca-pca-output := $(PRJ_DIR)/merged.eigenstrat.pca
smartpca-plot-output := $(PRJ_DIR)/merged.eigenstrat.plot
smartpca-eval-output := $(PRJ_DIR)/merged.eigenstrat.eval
smartpca-log-output := $(PRJ_DIR)/merged.eigenstrat.log

all: check-env create-project-dir $(smartpca-pca-output)

$(smartpca-pca-output) $(smartpca-plot-output) $(smartpca-eval-output) $(smartpca-log-output): $(genotype-output) $(snp-output) $(indiv-output)
	$(PERL) $(SMARTPCA) \
		-i $(genotype-output) \
		-a $(snp-output) \
		-b $(indiv-output) \
		-o $(smartpca-pca-output) \
		-p $(smartpca-plot-output) \
		-e $(smartpca-eval-output) \
		-l $(smartpca-log-output) \
		-k 10 \
		-m 0

$(genotype-output) $(snp-output) $(indiv-output): $(merge-par)
	$(CONVERTF) -p $(merge-par)

$(merge-par): $(INPUT_PED) $(INPUT_MAP)
	echo "genotypename: $(INPUT_PED)" > $(merge-par)
	echo "snpname:  $(INPUT_MAP)" >> $(merge-par)
	echo "indivname: $(INPUT_PED)" >> $(merge-par)
	echo "outputformat: EIGENSTRAT" >> $(merge-par)
	echo "genotypeoutname: $(genotype-output)" >> $(merge-par)
	echo "snpoutname: $(snp-output)" >> $(merge-par)
	echo "indivoutname: $(indiv-output)" >> $(merge-par)
	echo "familynames: NO" >> $(merge-par)

create-project-dir:
	mkdir -p $(PRJ_DIR)

check-env:
ifndef PRJ_DIR
	$(error PRJ_DIR is undefined!)
endif

ifndef INPUT_PED
	$(error INPUT_PED is undefined!)
endif

ifndef INPUT_MAP
	$(error INPUT_MAP is undefined!)
endif
