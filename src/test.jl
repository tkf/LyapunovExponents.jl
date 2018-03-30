module TestTools

using Base: rtoldefault
using Base.Test: @test, record, get_testset, Pass, Fail, Error, Broken
using DiffEqBase: DEProblem, set_u!, step!
using LyapunovExponents: LEProblem, dimension, get_integrator, current_state,
    PhaseTangentDynamics


macro test_nothrow(ex)
    quote
        @test begin
            $(esc(ex))
            true
        end
    end
end


function macro_kwargs(kwargs)
    params = map(
        kw -> begin
            @assert kw isa Expr
            @assert kw.head == :(=)
            Expr(:kw, kw.args...)
        end,
        kwargs)
    return Expr(:parameters, params...)
end


macro test_isapprox_elemwise(args...)
    quote
        record(get_testset(),
               @_test_isapprox_elemwise_result($(esc.(args)...)))
    end
end


macro _test_isapprox_elemwise_result(x, y, kwargs...)
    args = [
        esc(x), esc(y),
        QuoteNode(x), QuoteNode(y), QuoteNode(kwargs),
    ]
    Expr(:call,
         test_isapprox_elemwise,
         esc(macro_kwargs(kwargs)),
         args...)
end


function test_isapprox_elemwise(x, y, xexpr, yexpr, orig_kwargs;
                                skip::Bool = false,
                                broken::Bool = false,
                                io = STDOUT,
                                rtol::Real = rtoldefault(eltype(x), eltype(y)),
                                atol::Real = 0)
    appx = isapprox.(x, y; rtol=rtol, atol=atol)
    ok = all(appx)

    if !ok
        abserr = abs.(x .- y)
        amp = max.(abs.(x), abs.(y))
        th = atol .+ rtol .* amp
        relerr = abserr ./ amp

        print_with_color(:red, io, "Not pairwise isapprox")
        println(io, "  ", "rtol=", rtol, "  ", "atol=", atol)
        println(io, "x = ", xexpr)
        println(io, "y = ", yexpr)

        columns = [
            :i => 1:length(x),
            :x => x,
            :y => y,
            :relerr => relerr,
            :abserr => abserr,
            :th => th,
        ]
        table = []
        for (name, value) in columns
            push!(table, (name, sprint.(showcompact, value)))
        end
        mark = Dict(true => "✔", false => "!!")
        push!(table, ("ok?", [mark[x] for x in appx]))
        widths = [max(length(string(name)), maximum(length.(c)))
                  for (name, c) in table]

        for ((name, _), width) in zip(table, widths)
            print(io, "  ")
            print(io, lpad(name, width))
        end
        println(io)
        for i in 1:length(x)
            for ((_, c), width) in zip(table, widths)
                print(io, "  ")
                print(io, lpad(c[i], width))
            end
            println(io)
        end
    end

    orig_expr = Expr(:macrocall,
                     (:@test_isapprox_elemwise).args[1],
                     xexpr, yexpr, orig_kwargs...)

    return if skip
        Broken(:skipped, orig_expr)
    elseif broken
        if ok
            Error(:test_unbroken, orig_expr, ok, nothing)
        else
            Broken(:test, orig_expr)
        end
    else
        if ok
            Pass(:test, orig_expr, nothing, ok)
        else
            Fail(:test, orig_expr, nothing, ok)
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
                            t_evolve::Real = one(p1.tspan[1]),
                            kwargs...) where {P1 <: DEProblem,
                                              P2 <: DEProblem}

    for u in u_list
        i1 = get_integrator(p1)
        i2 = get_integrator(p2)
        set_u!(i1, copy(u))
        set_u!(i2, copy(u))
        step!(i1, t_evolve, true)
        step!(i2, t_evolve, true)
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
                            t_list = p1.tspan[1]:p1.tspan[1] + 1,
                            t_evolve::Real = one(p1.tspan[1]),
                            kwargs...) where {P1 <: DEProblem,
                                              P2 <: DEProblem}
    u_list = gen_u_list(rng, u_gen, num_u, p1.u0)
    if evolve
        test_same_dynamics(p1, p2, u_list;
                           t_evolve = t_evolve,
                           kwargs...)
    else
        test_same_dynamics(p1.f, p2.f, u_list, p_list, t_list; kwargs...)
    end
end


function test_tangent_dynamics_against_autodiff(
        prob::LEProblem, args...;
        kwargs...)
    @assert ! (prob.tangent_prob.f isa PhaseTangentDynamics)
    prob_ad = LEProblem(prob.phase_prob;
                        # t_attr is unused in the following, but it's
                        # a required argument:
                        t_attr = prob.t_attr,
                        # Q0 is required to match the size of
                        # tangent_prob.u0:
                        Q0 = prob.tangent_prob.u0[:, 2:end])
    prob_ad.tangent_prob.u0 .= prob.tangent_prob.u0
    prob_ad.tangent_prob.f :: PhaseTangentDynamics
    test_same_dynamics(prob.tangent_prob,
                       prob_ad.tangent_prob,
                       args...;
                       t_evolve = prob.t_renorm,
                       t_list = 0:prob.t_renorm,
                       kwargs...)
end

end  # module TestTools

using .TestTools
