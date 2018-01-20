using Base.Test
using LyapunovExponents: EXAMPLES, solve, lyapunov_exponents

@time @testset "Example $(ex.name)" for ex in [f() for f in EXAMPLES]
    solver = solve(ex; dim_lyap=length(ex.known_exponents))
    @test isapprox(ex.known_exponents, lyapunov_exponents(solver);
                   rtol=5e-2)   # TODO: Improve!
end
