using Base.Test
using LyapunovExponents: EXAMPLES, dimension, solve, lyapunov_exponents,
    LEProblem
using LyapunovExponents.Test: test_tangent_dynamics_against_autodiff

@time @testset "Tangent dynamics $(ex.name)" for ex in
        [ex for ex in [f() for f in EXAMPLES]
         if ex.tangent_dynamics! != nothing]
    num_u = 5  # number of random states at which the systems are tested
    test_tangent_dynamics_against_autodiff(LEProblem(ex), num_u)
end

@time @testset "Example $(ex.name)" for ex in [f() for f in EXAMPLES]
    for dim_lyap in 1:dimension(ex)
        if dim_lyap != length(ex.known_exponents) && (
                contains(ex.name, "Linz & Sprott (1999)") ||
                contains(ex.name, "standard map"))
            # TODO: check why they don't work
            continue
        end
        solver = solve(ex; dim_lyap=dim_lyap)
        @show dim = min(dim_lyap, length(ex.known_exponents))
        @test isapprox((@show ex.known_exponents[1:dim]),
                       (@show lyapunov_exponents(solver)[1:dim]);
                       rtol=ex.rtol, atol=ex.atol)
    end
end
