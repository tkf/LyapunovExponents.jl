using PkgBenchmark: judge, target_result, baseline_result,
    writeresults, export_markdown

mkpath("tmp")
JUDGE_BASELINE = get(ENV, "JUDGE_BASELINE", "^HEAD")
results = judge("LyapunovExponents", JUDGE_BASELINE)
writeresults("tmp/benchmarks_target.json", target_result(results))
writeresults("tmp/benchmarks_baseline.json", baseline_result(results))
export_markdown("tmp/benchmarks_judge.md", results)
showall(results)
println()
