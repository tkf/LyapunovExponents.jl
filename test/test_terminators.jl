using Base.Test
using LyapunovExponents
using LyapunovExponents: report, dimension, arnold_cat_map, van_der_pol

@time @testset "Terminate stable system $(ex.name)" for ex in [
        arnold_cat_map(M = [0.5 0.1; 1 0.1],
                       t_attr = 100000).example,
        van_der_pol(params = [:a => 0.1],
                    t_attr = 100000).example,
        ]
    for dim_lyap in 1:dimension(ex)
        println("$(ex.name) dim_lyap=$dim_lyap")
        @time sol = solve(ex; record=true, dim_lyap=dim_lyap)
        report(sol)
        @test lyapunov_exponents(sol)[1] < 0
        @test sol.converged
    end
end
