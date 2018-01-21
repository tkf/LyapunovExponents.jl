JULIA = julia --color=yes
DOCS_PORT = 35729
LIVERELOAD = livereload --port $(DOCS_PORT) --wait 0.5

.PHONY: *

test:
	$(JULIA) --check-bounds=yes test/runtests.jl

docs:
	$(JULIA) docs/make.jl
	cd docs && mkdocs build

docs-cont:
	rg --files | entr $(MAKE) docs

serve:
	cd docs/site && $(LIVERELOAD)
#	cd docs && mkdocs serve --dev-addr=localhost:$(DOCS_PORT)
# Somehow, "mkdocs serve" version didn't work well; it didn't update
# the page on-the-fly.  That's why I'm using livereload here.

docs-serve:
	$(MAKE) docs-cont &
	$(MAKE) serve
