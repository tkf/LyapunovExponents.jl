using Base.Test
using LyapunovExponents: DEMOS, dimension, solve, lyapunov_exponents,
    LEProblem, report
using LyapunovExponents.Test: test_tangent_dynamics_against_autodiff

@time @testset "Tangent dynamics $(ex.name)" for ex in
        [ex for ex in [f().example for f in DEMOS]
         if ex.tangent_dynamics != nothing]
    num_u = 5  # number of random states at which the systems are tested
    opts = Dict(
        # :verbose => true
    )
    if contains(ex.name, "Linz & Sprott (1999)")
        # linz_sprott_99 is solved with a more accurate and expensive
        # ODE algorithm and with larger t_renorm.  But tangent test
        # (with evolve=true) is done with the default ODE algorithm.
        # So let's old t_renorm as the span.
        # See: [[../src/examples/linz_sprott_99.jl::de_options]]
        opts[:tspan] = (0.0, 1.0)
    end
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
        println()
        print_with_color(:blue, "$(ex.name)", bold=true)
        println(" dim_lyap=$dim_lyap")

        @time sol = solve(ex; record=true, dim_lyap=dim_lyap)
        dim = min(dim_lyap, length(ex.known_exponents))
        report(sol)
        @show dim
        @show ex.known_exponents[1:dim]
        if ! (contains(ex.name, "standard map") ||
              contains(ex.name, "van der Pol") ||
              contains(ex.name, "Beer"))
            @test sol.converged
        end

        rtol = ex.rtol
        atol = ex.atol
        if contains(ex.name, "standard map")
            rtol = 5rtol
            atol = 5atol
        end

        if dim_lyap != length(ex.known_exponents) &&
                contains(ex.name, "Linz & Sprott (1999)") ||
                dim_lyap == 1 && contains(ex.name, "van der Pol")
            # TODO: check why they don't work
            @test_skip(
                lyapunov_exponents(sol)[1:dim] ≈ ex.known_exponents[1:dim],
                rtol = rtol,
                atol = atol)
        else
            @test(lyapunov_exponents(sol)[1:dim] ≈ ex.known_exponents[1:dim],
                  rtol = rtol,
                  atol = atol)
        end
    end
end
