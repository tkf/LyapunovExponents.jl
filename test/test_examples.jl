using Base.Test
using LyapunovExponents: EXAMPLES, solve, lyapunov_exponents

@time @testset "Example $(ex.name)" for ex in [f() for f in EXAMPLES]
    solver = solve(ex; dim_lyap=length(ex.known_exponents))
    @test isapprox((@show ex.known_exponents),
                   (@show lyapunov_exponents(solver));
                   rtol=ex.rtol, atol=ex.atol)
end
