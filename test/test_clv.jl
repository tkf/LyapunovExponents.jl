using Base.Test
using IterTools: product
using LyapunovExponents
using LyapunovExponents: DEMOS, objname, dimension
using LyapunovExponents.CovariantVectors: goto!

@time @testset "CLV: $(objname(f))" for f in DEMOS
    @testset "num_backward_tran=$nbt brx=$test_brx dim_lyap=$dim_lyap" for (
                nbt, test_brx, dim_lyap,
            ) in product(
                0:2,            # num_backward_tran
                [true, false],  # test_brx
                2:dimension(f().example),  # dim_lyap
            )
        prob = CLVProblem(f(dim_lyap=dim_lyap).prob::LEProblem;
                          num_clv = 5,
                          num_forward_tran = nbt * 11,  # for extra variation
                          num_backward_tran = nbt)

        dims = (dp, dl) = size(prob.Q0)
        @assert dims == (dimension(prob.phase_prob), dim_lyap)
        solver = init(prob;
                      backward_dynamics = CLV.BackwardDynamicsWithD)

        forward = forward_dynamics!(solver)
        @test length(forward) == prob.num_clv + prob.num_backward_tran
        R_prev = [Matrix{Float64}(dl, dl) for _ in 1:length(forward)]
        G = [Matrix{Float64}(dp, dl) for _ in 1:length(forward)]
        M = [Matrix{Float64}(dp, dp) for _ in 1:length(forward)]

        # Recall: 𝑮ₙ₊ₖ 𝑪ₙ₊ₖ 𝑫ₖ,ₙ = 𝑴ₖ,ₙ 𝑮ₙ 𝑪ₙ = 𝑮ₙ₊ₖ 𝑹ₖ,ₙ 𝑪ₙ
        # (Eq. 32, Ginelli et al., 2013)

        for (n, _) in indexed_forward_dynamics!(forward)
            G[n] .= CLV.G(forward)  # 𝑮ₙ
            M[n] .= CLV.M(forward)  # 𝑴ₖ,ₙ
            R_prev[n] .= CLV.R_prev(forward)  # 𝑹ₖ,ₙ
        end
        @testset "forward.R_history[$n]" for n in 1:length(forward.R_history)
            @test forward.R_history[n] == R_prev[n]
        end
        @assert forward.R_history == R_prev

        if test_brx
            brx = goto!(solver, CLV.BackwardRelaxer)
            @testset "brx.R[$n]" for n in 1:length(brx.R_history)
                @test brx.R_history[n] == R_prev[n]
            end
            @assert all(brx.R_history .== R_prev)
        end

        backward = backward_dynamics!(solver)
        num_clv = length(backward)
        @assert backward.R_history == R_prev[1:num_clv]
        C = [Matrix{Float64}(dl, dl) for _ in 1:num_clv]
        D = [Matrix{Float64}(dl, dl) for _ in 1:num_clv]
        C[end] .= CLV.C(backward)
        for (n, Cn) in indexed_backward_dynamics!(backward)
            @test CLV.R(backward) == R_prev[n+1]  # 𝑹ₖ,ₙ
            C[n] .= Cn               # 𝑪ₙ
            D[n] .= CLV.D(backward)  # 𝑫ₖ,ₙ
        end

        # @testset "𝑮ₙ₊ₖ 𝑹ₖ,ₙ ≈ 𝑴ₖ,ₙ 𝑮ₙ (n=$n)" for n in 1:num_clv-1
        @testset "Gₙ₊ₖ Rₖ,ₙ ≈ Mₖ,ₙ Gₙ (n=$n)" for n in 1:num_clv-1
            # TODO: improve rtol
            @test G[n+1] * R_prev[n+1] ≈ M[n] * G[n]  rtol=5e-2
            # v-- commutative diagram (1)
        end
        # @testset "𝑪ₙ₊ₖ 𝑫ₖ,ₙ ≈ 𝑹ₖ,ₙ 𝑪ₙ (n=$n)" for n in 1:num_clv-1
        @testset "Cₙ₊ₖ Dₖ,ₙ ≈ Rₖ,ₙ Cₙ (n=$n)" for n in 1:num_clv-1
            Rn = R_prev[n+1]
            @test C[n+1] * D[n] ≈ Rn * C[n]
            # v-- commutative diagram (2)
        end
        # @testset "𝑮ₙ₊ₖ 𝑪ₙ₊ₖ 𝑫ₖ,ₙ = 𝑴ₖ,ₙ 𝑮ₙ 𝑪ₙ (n=$n)" for n in 1:num_clv-1
        @testset "Gₙ₊ₖ Cₙ₊ₖ Dₖ,ₙ = Mₖ,ₙ Gₙ Cₙ (n=$n)" for n in 1:num_clv-1
            # TODO: improve rtol
            @test G[n+1] * C[n+1] * D[n] ≈ M[n] * G[n] * C[n]  rtol=1e-1
        end

        #  ─── M[n] ──▶
        # ▲                 ▲
        # │                 │
        # G[n]    (1)      G[n+1]
        # │                 │
        #  ─── R[n] ──▶
        # ▲                 ▲
        # │                 │
        # C[n]    (2)      C[n+1]
        # │                 │
        #  ─── D[n] ──▶

    end
end
