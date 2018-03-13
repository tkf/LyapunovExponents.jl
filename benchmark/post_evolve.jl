using PkgBenchmark
using DiffEqBase: DiscreteProblem
using LyapunovExponents: TangentRenormalizer, MLERenormalizer,
    post_evolve!, get_integrator, objname

""" Fake LESolver for benchmarking `post_evolve!`. """
function noop_tangent_integrator(N, dim_lyap)
    tangent_prob = DiscreteProblem((x...) -> nothing,  # irrelevant
                                   zeros(N, 1+ dim_lyap),
                                   (0, 1))  # irrelevant
    return get_integrator(tangent_prob)
end

suite = BenchmarkGroup()
for renormalizer in [TangentRenormalizer, MLERenormalizer]
    bg = suite[objname(renormalizer)] = BenchmarkGroup()
    for N in 2.^(6:7), dim_lyap in (renormalizer === TangentRenormalizer ?
                                    floor.(Int, [0.1, 0.5, 0.9, 1] .* N) :
                                    [1])
        stage = renormalizer(noop_tangent_integrator(N, dim_lyap),
                             1, 1)
        bg["N=$N", "dim_lyap=$dim_lyap"] =
            @benchmarkable(post_evolve!($stage),
                           setup = (randn!($stage.tangent_state)))
    end
end

suite
