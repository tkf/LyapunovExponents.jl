using Base.Test
using LyapunovExponents: DEMOS, dimension, solve, lyapunov_exponents,
    LEProblem
using LyapunovExponents.Test: test_tangent_dynamics_against_autodiff

@time @testset "Tangent dynamics $(ex.name)" for ex in
        [ex for ex in [f().example for f in DEMOS]
         if ex.tangent_dynamics != nothing]
    num_u = 5  # number of random states at which the systems are tested
    opts = Dict(
        # :verbose => true
    )
    test_tangent_dynamics_against_autodiff(LEProblem(ex), num_u; opts...)
    if contains(ex.name, "Hénon map")
        # TODO: It fails maybe because some randomly generated states
        # are outside the domain of Hénon map.
        continue
    else
        test_tangent_dynamics_against_autodiff(LEProblem(ex), num_u;
                                               evolve = true,
                                               opts...)
    end
end

@time @testset "Example $(ex.name)" for ex in [f().example for f in DEMOS]
    for dim_lyap in 1:dimension(ex)
        println("$(ex.name) dim_lyap=$dim_lyap")
        @time sol = solve(ex; dim_lyap=dim_lyap)
        @show dim = min(dim_lyap, length(ex.known_exponents))
        @show ex.known_exponents[1:dim]
        @show lyapunov_exponents(sol)[1:dim]
        if dim_lyap != length(ex.known_exponents) &&
                contains(ex.name, "Linz & Sprott (1999)") ||
                dim_lyap == 1 && contains(ex.name, "van der Pol")
            # TODO: check why they don't work
            @test_skip isapprox(ex.known_exponents[1:dim],
                                lyapunov_exponents(sol)[1:dim];
                                rtol=ex.rtol, atol=ex.atol)
        else
            @test isapprox(ex.known_exponents[1:dim],
                           lyapunov_exponents(sol)[1:dim];
                           rtol=ex.rtol, atol=ex.atol)
        end
    end
end
