using Base.Test
using LyapunovExponents

@testset "User interface" begin
    args = [
        (x...) -> nothing,  # f
        [0.0],              # u0
        (0.0, 1.0),         # tspan
        nothing,            # p
    ]
    for lep in [DiscreteLEProblem, ContinuousLEProblem]
        @test begin
            lep(args...; num_attr=1)
            true
        end
        @test begin
            lep(args..., 1)
            true
        end
        try
            lep(args...)
            println("Calling $lep w/o num_attr didn't throw.")
            @test false
        catch err
            @test isa(err, ErrorException)
            @test contains(err.msg, "`num_attr` is required")
        end
    end
end
