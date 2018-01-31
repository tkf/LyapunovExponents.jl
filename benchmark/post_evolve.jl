using PkgBenchmark
using DiffEqBase: DiscreteProblem
using LyapunovExponents: LESolver, post_evolve!, get_integrator

""" Fake LESolver for benchmarking `post_evolve!`. """
function fake_le_solver(N, dim_lyap)
    tangent_prob = DiscreteProblem((x...) -> nothing,  # irrelevant
                                   zeros(N, 1+ dim_lyap),
                                   (0, 1))  # irrelevant
    return LESolver(get_integrator(tangent_prob))
end

@benchgroup "post_evolve!(::LESolver)" begin
    for N in 2.^(6:7), dim_lyap in floor.(Int, [0.1, 0.5, 0.9, 1] .* N)
        solver = fake_le_solver(N, dim_lyap)
        @bench(("N=$N", "dim_lyap=$dim_lyap"),
               post_evolve!($solver),
               setup=(randn!($solver.tangent_state)))
    end
end
