using Base.Test
using OnlineStats: Series, fit!
using LyapunovExponents: VecMean, VecVariance

@time @testset "OnlineStats" begin
    max_dim = 10
    len = 5
    rng = MersenneTwister(0)
    @testset "VecMean dim=$dim" for dim in 1:max_dim
        stat = VecMean(dim)
        os = Series(stat)
        xs = randn(rng, dim, len)
        for i in 1:len
            fit!(os, xs[:, i])
            @test mean(stat) ≈ mean(xs[:, 1:i], 2)[:, 1]
        end
    end
    @testset "VecVariance dim=$dim" for dim in 1:max_dim
        stat = VecVariance(dim)
        os = Series(stat)
        xs = randn(rng, dim, len)
        for i in 1:len
            fit!(os, xs[:, i])
            @test mean(stat) ≈ mean(xs[:, 1:i], 2)[:, 1]
            @test var(stat) ≈ var(xs[:, 1:i], 2)[:, 1] nans=true
        end
    end
end
