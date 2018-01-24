using Base.Test
using DifferentialEquations: ODEProblem, solve
using LyapunovExponents
using LyapunovExponents: get_integrator, keepgoing!, LEProblem

abserr(x, y) = @views maximum(abs.(x[:] .- y[:]))

""" Symmetric mean absolute relative error (like sMAPE but w/o 100x) """
smarerr(x, y) = 2 * @views maximum(abs.(x[:] .- y[:]) ./ (x[:] .+ y[:]))

@time @testset "Keepgoing $(ex.name)" for ex in [
        LyapunovExponents.lorenz_63().example,
        LyapunovExponents.linz_sprott_99().example,
        LyapunovExponents.van_der_pol().example,
        ]
    ode = LEProblem(ex).phase_prob

    points = 100
    t0, t1 = ode.tspan
    t2 = (t1 - t0) * 2 + t0

    long_ode = ODEProblem(ode.f, copy(ode.u0), (t0, t2), ode.p)
    long_sol = solve(long_ode)
    times = collect(linspace(t0, t2, points))
    long_data = hcat(map(long_sol, times)...)

    integrator = get_integrator(ode; save_everystep=true)
    keepgoing!(integrator)
    first_data = [integrator.sol(t) for t in times if t < t1]
    half_sol_in_first = copy(integrator.sol[end])

    keepgoing!(integrator)
    second_data = [integrator.sol(t - t1) for t in times if t >= t1]
    half_sol_in_second = integrator.sol(t0)

    cat_data = hcat(vcat(first_data, second_data)...)

    @show abserr(cat_data, long_data)
    @show smarerr(cat_data, long_data)
    @test isapprox(cat_data, long_data;
                   rtol=integrator.opts.reltol,
                   atol=integrator.opts.abstol)

    @test half_sol_in_first ≈ half_sol_in_second
end
