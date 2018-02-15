JULIA = julia --color=yes
DOCS_PORT = 35729
LIVERELOAD = livereload --port $(DOCS_PORT) --wait 0.5

.PHONY: *

test:
	$(JULIA) --check-bounds=yes test/runtests.jl

plot replot: %:
	$(MAKE) --directory=docs $*

judge:
	OMP_NUM_THREADS=1 $(JULIA) -e 'using PkgBenchmark: judge; showall(judge("LyapunovExponents", "$(JUDGE_BASELINE)"; promptsave=false, promptoverwrite=false)); println()'
JUDGE_BASELINE = HEAD^

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
