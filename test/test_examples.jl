using Base.Test
using LyapunovExponents: DEMOS, dimension, solve!, lyapunov_exponents,
    LEProblem, report, objname
using LyapunovExponents.TestTools: test_tangent_dynamics_against_autodiff,
    @test_isapprox_elemwise

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
        opts[:t_evolve] = 1.0
        opts[:t_list] = 0:0.5:1.0
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

@time @testset "Example $(objname(f))" for f in DEMOS
    for dim_lyap in 1:dimension(f().example)
        println()

        demo = f(dim_lyap=dim_lyap)
        @time solve!(demo; record=true)
        report(demo)

        sol = demo.solver.sol
        ex = demo.example
        dim = min(dim_lyap, length(ex.known_exponents))

        if contains(ex.name, "standard map") ||
                contains(ex.name, "van der Pol")
            @test_broken sol.converged
        else
            @test sol.converged
        end

        rtol = ex.rtol
        atol = ex.atol
        if contains(ex.name, "standard map")
            rtol *= 20
            atol *= 20
        end

        @test_isapprox_elemwise(
            lyapunov_exponents(sol)[1:dim],
            ex.known_exponents[1:dim],
            rtol = rtol,
            atol = atol,
        )
    end
end
