using Documenter, LyapunovExponents

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-alabaster"),
    repo   = "github.com/tkf/LyapunovExponents.jl.git",
    julia  = "0.6",
)
