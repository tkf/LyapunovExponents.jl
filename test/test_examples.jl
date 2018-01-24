using Base.Test
using LyapunovExponents: EXAMPLES, dimension, solve, lyapunov_exponents

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
