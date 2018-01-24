using Base.Test
using LyapunovExponents: EXAMPLES, solve, lyapunov_exponents

@time @testset "Example $(ex.name)" for ex in [f() for f in EXAMPLES]
    for dim_lyap in [length(ex.known_exponents), 1]
        if dim_lyap == 1 && (
                contains(ex.name, "Linz & Sprott (1999)") ||
                contains(ex.name, "standard map"))
            # TODO: check why they don't work with MLESolver
            continue
        end
        solver = solve(ex; dim_lyap=dim_lyap)
        @test isapprox((@show ex.known_exponents[1:dim_lyap]),
                       (@show lyapunov_exponents(solver));
                       rtol=ex.rtol, atol=ex.atol)
    end
end
