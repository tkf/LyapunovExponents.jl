using Documenter, LyapunovExponents

using Plots
gr()

function should_plot(env = ENV)
    (@show get(env, "TRAVIS_BRANCH", "")) == "master" &&
        (@show get(env, "TRAVIS_PULL_REQUEST", "false")) == "false"
end

# Generate src/gallery/examples/*.{md,png}.  In principle, I can
# generate those files in build/gallery/examples/ after `makedocs()`
# below.  However, resolving @ref etc. are done by `makedocs()` so
# it's better to generate Markdown files here.
cd(dirname(@__FILE__)) do
    run(`make gallery_md`)
    if should_plot()
        include("plot_gallery.jl")
    end
end

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-alabaster"),
    repo   = "github.com/tkf/LyapunovExponents.jl.git",
    julia  = "0.6",
)
