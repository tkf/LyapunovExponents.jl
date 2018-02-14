module Test

using Base: rtoldefault
using Base.Test: @test
using DiffEqBase: DEProblem
using LyapunovExponents: LEProblem, dimension, phase_tangent_state,
    get_tangent_prob, get_integrator, de_prob, keepgoing!, current_state


macro test_nothrow(ex)
    quote
        @test begin
            $(esc(ex))
            true
        end
    end
end


macro display(ex)
    quote
        println($("$ex ="))
        display($(esc(ex)))
        println()
    end
end


function gen_u_list(rng, u_gen, num_u::Int, u::T) where {T}
    u_list = T[]
    for _ in 1:num_u
        push!(u_list, u_gen(rng, eltype(u), size(u)))
    end
    return u_list
end


function rtol_exponents_default(x, y, n=50)
    return linspace(log10(rtoldefault(eltype(x), eltype(y))), 0.5, n)
end


function analyze_diff(x, y, ratio=0.1;
                      rtol_exponents = rtol_exponents_default(x, y))
    for n in rtol_exponents
        rtol = 10^n
        bad = mean(@. abs(x - y) > rtol * max(abs(x), abs(y))) / length(x)
        if bad < ratio
            return rtol, bad
        end
    end
    0, 0
end


macro show_diff_info(x, y)
    quote
        rtol, bad = analyze_diff($(esc(x)), $(esc(y)))
        if bad > 0
            @printf("%.1f%% elements are %s ≉ %s with relative tolerance %g\n",
                    bad * 100, $(string(x)), $(string(y)), rtol)
        end
    end
end


function compare_states(u1, u2; verbose=false, kwargs...)
    if verbose
        @display u1
        @display u2
        @display u1 - u2
    end
    @show_diff_info(u1, u2)
    @test isapprox(u1, u2; kwargs...)
end


function test_same_dynamics(f1!, f2!, u_list, p_list, t_list; kwargs...)
    for u in u_list
        for p in p_list
            for t in t_list
                u1 = similar(u)
                u2 = similar(u)
                f1!(u1, copy(u), p, t)
                f2!(u2, copy(u), p, t)
                compare_states(u1, u2; kwargs...)
            end
        end
    end
end


function test_same_dynamics(p1::P1, p2::P2, u_list::AbstractArray;
                            kwargs...) where {P1 <: DEProblem,
                                              P2 <: DEProblem}
    @assert p1.tspan == p2.tspan

    i1 = get_integrator(p1)
    i2 = get_integrator(p2)
    for u in u_list
        keepgoing!(i1, copy(u))
        keepgoing!(i2, copy(u))
        u1 = current_state(i1)
        u2 = current_state(i2)
        compare_states(u1, u2; kwargs...)
    end
end


function test_same_dynamics(p1::P1, p2::P2, num_u::Int;
                            evolve = false,
                            rng = Base.GLOBAL_RNG,
                            u_gen = randn,
                            p_list = [p1.p],
                            kwargs...) where {P1 <: DEProblem,
                                              P2 <: DEProblem}
    u_list = gen_u_list(rng, u_gen, num_u, p1.u0)
    if evolve
        test_same_dynamics(p1, p2, u_list; kwargs...)
    else
        t_list = p1.tspan[1]:p1.tspan[2]
        test_same_dynamics(p1.f, p2.f, u_list, p_list, t_list; kwargs...)
    end
end


function test_tangent_dynamics_against_autodiff(
        prob::LEProblem, args...;
        dim_lyap = dimension(prob.phase_prob),
        kwargs...)
    @assert prob.tangent_dynamics! != nothing
    prob_ad = LEProblem(prob.phase_prob, prob.num_attr)
    u0 = phase_tangent_state(prob)
    test_same_dynamics(get_tangent_prob(prob, u0),
                       get_tangent_prob(prob_ad, u0),
                       args...; kwargs...)
end

end

using .Test
