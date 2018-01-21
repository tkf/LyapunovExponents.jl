JULIA = julia --color=yes

.PHONY: *

test:
	$(JULIA) --check-bounds=yes test/runtests.jl

docs:
	$(JULIA) docs/make.jl
