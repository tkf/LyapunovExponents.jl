using Base.Test
using LyapunovExponents: is_semi_unitary

@testset "Utils: is_semi_unitary" begin
    @test is_semi_unitary(eye(3))
    @test ! is_semi_unitary(zeros(Int, 2, 2))
    @test ! is_semi_unitary(ones(Int, 2, 2))
    @test is_semi_unitary(qr([[1, 2] [2, 1]])[1])
end
