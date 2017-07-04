.PHONY: init test clean docs check-sphinx

SHELL := /bin/bash
SPHINX := $(shell which sphinx-build)
GH_PAGES_SOURCES = docs/Makefile docs/source
MASTER_HEAD_COMMIT = $(shell git log master -1 --pretty=short --abbrev-commit)

init:
	pip install -r requirements.txt

test:
	echo "Please implement me!"

clean:
	find ./yaps2 -name "*.pyc" -exec rm {} \;

docs: check-sphinx
	git checkout gh-pages
	rm -rf build _sources _static
	git checkout doc-setup $(GH_PAGES_SOURCES)
	git reset HEAD
	$(MAKE) -f docs/Makefile html
	mv -fv build/html/* ./
	rm -rf $(GH_PAGES_SOURCES) build
	git add -A
	git ci -m \
		"Generated gh-pages for $(MASTER_HEAD_COMMIT)"

check-sphinx:
ifndef SPHINX
	$(error "sphinx not installed:  run 'pip install Sphinx'")
endif
