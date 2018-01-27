using Base.Test
using Plots
using LyapunovExponents

@time @testset "Smoke test demo" begin
    @test begin
        demo = solve!(LyapunovExponents.lorenz_63(num_attr=3))
        plot(demo, show = false)
        true
    end
end
