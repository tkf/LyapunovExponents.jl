using PkgBenchmark: benchmarkpkg, writeresults, export_markdown

mkpath("tmp")
results = benchmarkpkg("LyapunovExponents")
writeresults("tmp/benchmarks_current.json", results)
export_markdown("tmp/benchmarks_current.md", results)
showall(results)
println()
