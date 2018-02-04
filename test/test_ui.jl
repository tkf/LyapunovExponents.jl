using Base.Test
using LyapunovExponents
using LyapunovExponents.Test: @test_nothrow

@testset "User interface" begin
    args = [
        (x...) -> nothing,  # f
        [0.0],              # u0
        (0.0, 1.0),         # tspan
        nothing,            # p
    ]
    for lep in [DiscreteLEProblem, ContinuousLEProblem]
        @test_nothrow lep(args...; num_attr=1)
        @test_nothrow lep(args..., 1)
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
