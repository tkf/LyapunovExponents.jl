using Base.Test
using LyapunovExponents
using LyapunovExponents.Test: @test_nothrow

@testset "User interface" begin
    args = [
        (x...) -> nothing,  # f
        [0.0],              # u0
        1,                  # t_renorm
        nothing,            # p
    ]
    for lep in [DiscreteLEProblem, ContinuousLEProblem]
        @test_nothrow lep(args...; t_attr=1)
        @test_nothrow lep(args..., 1)
        try
            lep(args...)
            println("Calling $lep w/o t_attr didn't throw.")
            @test false
        catch err
            @test isa(err, ErrorException)
            @test contains(err.msg, "`t_attr` is required")
        end
    end
end
