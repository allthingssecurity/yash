PYTHON ?= python

.PHONY: help verify reproduce test
help:
	@echo 'Start: README.md and docs/usage.md'
	@echo 'make verify: check recorded lead Design #3 results'
	@echo 'make reproduce: replay lead Design #3 from recorded evidence, locally'
	@echo 'Make another candidate: yash-vaccine init NAME'
	@echo 'Inspect commands: yash-vaccine --help'

verify:
	$(PYTHON) -m yash_vaccine verify

reproduce:
	$(PYTHON) -m yash_vaccine reproduce

test:
	$(PYTHON) -m unittest discover -s tests -v
