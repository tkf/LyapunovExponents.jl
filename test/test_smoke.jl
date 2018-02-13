using Base.Test
using Plots
using OnlineStats: CovMatrix
using LyapunovExponents
using LyapunovExponents: VecVariance
using LyapunovExponents.Test: @test_nothrow

@time @testset "Smoke test CLV" begin
    @test_nothrow begin
        prob = LyapunovExponents.lorenz_63(num_attr=3).prob :: LEProblem
        solver = CLVSolver(prob)
        solve!(solver)
    end
    @testset "Recording CLV" begin
        prob = LyapunovExponents.lorenz_63(num_attr=3).prob :: LEProblem
        solver = CLVSolver(
            prob;
            record = [:G, :C, :D, :x])
        solve!(solver)
        @test isdefined(solver.sol, :G_history)
        @test isdefined(solver.sol, :C_history)
        @test isdefined(solver.sol, :D_history)
        @test isdefined(solver.sol, :x_history)
        num_fwd = solver.prob.num_clv + solver.prob.num_backward_tran
        @test length(solver.sol.G_history) == num_fwd
        @test length(solver.sol.C_history) == solver.prob.num_clv
        @test length(solver.sol.D_history) == solver.prob.num_clv
        @test length(solver.sol.x_history) == solver.prob.num_clv
    end
end

@time @testset "Smoke test demo" begin
    @test begin
        demo = solve!(LyapunovExponents.lorenz_63(num_attr=3))
        plot(demo, show = false)
        true
    end
    @test begin
        demo = solve!(LyapunovExponents.lorenz_63(num_attr=3, dim_lyap=1))
        plot(demo, show = false)
        true
    end
    @testset "main_stat = $ST" for ST in [VecVariance, CovMatrix]
        demo = LyapunovExponents.lorenz_63(num_attr=3)
        @test_nothrow solve!(demo; main_stat = ST)
        @test isa(demo.solver.solver.main_stat, ST)
        @test_nothrow plot(demo, show = false)
    end
end
