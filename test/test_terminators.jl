using Base.Test
using LyapunovExponents
using LyapunovExponents: stable_le_error, FixedPointConvDetail,
    NonNegativeAutoCovConvDetail, NonPeriodicConvDetail,
    report, dimension, arnold_cat_map, van_der_pol

@time @testset "stable_le_error" begin
    with_negative(xs) = [xs, -xs]

    for xs in with_negative(repeat(1:3, outer=10))
        err, detail = stable_le_error(xs)
        @test detail.period == 3
    end

    for xs in with_negative(repeat([1], outer=10))
        err, detail = stable_le_error(xs)
        @test detail isa FixedPointConvDetail
        @test detail.var == 0
        @test err == 0
    end

    for xs in with_negative(collect(1:20))
        err, detail = stable_le_error(xs, cutoff=5)
        @test detail isa NonNegativeAutoCovConvDetail
        @test err == Inf
    end

    for xs in with_negative(collect(1:20))
        err, detail = stable_le_error(xs)
        @test detail isa NonPeriodicConvDetail
        @test err == Inf
    end
end

@time @testset "Terminate stable system $(ex.name)" for ex in [
        arnold_cat_map(M = [0.5 0.1; 1 0.1],
                       t_attr = 100000).example,
        van_der_pol(params = [:a => 3],
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
