EXAMPLE_SCRIPTS = $(sort $(wildcard src/gallery/examples/*.jl))
EXAMPLE_MD = $(patsubst %.jl, %.md, $(EXAMPLE_SCRIPTS))

.PHONY: gallery_md

gallery_md: $(EXAMPLE_MD)

$(EXAMPLE_MD): %.md: %.jl gen_example_md
	./gen_example_md $< > $@

mkdocs.yml: $(EXAMPLE_MD)
	sed '/- Gallery:/q' $@ > mkdocs.yml.tmp
	bash -c "for name in $(patsubst build/%.md, %.md, $(EXAMPLE_MD)); do echo \"      - \$$name\"; done" >> mkdocs.yml.tmp
	mv mkdocs.yml.tmp mkdocs.yml
