using Documenter, LyapunovExponents

# Generate src/gallery/examples/*.{md,png}.  In principle, I can
# generate those files in build/gallery/examples/ after `makedocs()`
# below.  However, resolving @ref etc. are done by `makedocs()` so
# it's better to generate Markdown files here.
cd(dirname(@__FILE__)) do
    run(`make gallery_md`)
    include("plot_gallery.jl")
end

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-alabaster"),
    repo   = "github.com/tkf/LyapunovExponents.jl.git",
    julia  = "0.6",
)
