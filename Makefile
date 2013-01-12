# Makefile for IMUSim documentation and tests

.PHONY: default doc api tutorial test

export PYTHONPATH:=$(shell pwd):${PYTHONPATH}

default:
	@echo "This Makefile is only for building documentation and running tests."
	@echo "To build and install IMUSim, run:"
	@echo
	@echo "    python setup.py build"
	@echo "    python setup.py install (as root or administrator)"
	@echo
	@echo "Or, from the parent directory, with setuptools installed:"
	@echo
	@echo "    easy_install imusim"
	@echo
	@echo "This will also download and install required dependencies."

c:
	cython `find . -name \*.pyx`

doc: api sphinx

api:
	mkdir -p docs/build/html/api/
	epydoc -v --docformat epytext --no-frames --no-private \
		--graph classtree --exclude ^imusim.tests --exclude ^imusim.all \
		imusim -o docs/build/html/api/ -n IMUSim -u http://www.imusim.org/

sphinx:
	cd docs; sphinx-build -b html . build/html/

test:
	nosetests -v -x --with-coverage --cover-package=imusim --cover-inclusive \
		--cover-html --cover-html-dir docs/build/html/coverage --with-id \
		imusim/tests
