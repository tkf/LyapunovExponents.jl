EXAMPLE_SCRIPTS = $(sort $(wildcard src/gallery/examples/*.jl))
EXAMPLE_MD = $(patsubst %.jl, %.md, $(EXAMPLE_SCRIPTS))
GALLERY_MD = src/gallery/Ginelli_2007_Figure_1a.md \
src/gallery/Lorenz_CLV.md

.PHONY: gallery_md

gallery_md: $(EXAMPLE_MD) $(GALLERY_MD) src/gallery/index.md

$(EXAMPLE_MD): %.md: %.jl gen_example_md
	./gen_example_md $< > $@

$(GALLERY_MD): %.md: %.jl %.md.head gen_gallery_md
	./gen_gallery_md $*.md.head $*.jl > $@

src/gallery/index.md: $(GALLERY_MD) $(EXAMPLE_MD) gen_gallery_index_md
	./gen_gallery_index_md $(GALLERY_MD) $(EXAMPLE_MD) > $@

mkdocs.yml: src/gallery/index.md $(GALLERY_MD) $(EXAMPLE_MD)
	sed '/- Gallery:/q' $@ > mkdocs.yml.tmp
	bash -c "for name in $(patsubst src/%.md, %.md, $^); do echo \"      - \$$name\"; done" >> mkdocs.yml.tmp
	mv mkdocs.yml.tmp mkdocs.yml


JULIA = julia --color=yes

plot:
	$(JULIA) plot_gallery.jl

replot:
	find src/gallery -name *.png -print0 | xargs -0 --no-run-if-empty rm -v
	$(MAKE) plot
