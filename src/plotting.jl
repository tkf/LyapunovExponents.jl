using RecipesBase
using StatsBase: autocor


struct ConvErrorBars
    convergence::ConvergenceHistory
    exponents::AbstractVector
    i::Int
end

@recipe function f(data::ConvErrorBars;
                   threshold_color = :green,
                   unstable_color = :darkblue,
                   stable_color = :darkred)
    history = data.convergence
    exponents = data.exponents
    i = data.i

    kind_to_color = Dict(
        UnstableConvError => unstable_color,
        StableConvError => stable_color,
    )
    errorcolor = [kind_to_color[k] for k in history.kinds]
    thresholdcolor = [threshold_color for _ in errorcolor]
    invisibles = [nothing for _ in errorcolor]

    seriestype := :scatter
    markersize --> [3 5]
    label := ""
    markershape --> :+
    # markerstrokestyle := [:solid :dash]
    markercolor := [errorcolor invisibles]
    markerstrokecolor := [errorcolor thresholdcolor]
    yerror := [history.errors[i] history.thresholds[i]]
    x = history.orth
    l = [exponents[n] for n in history.orth]
    y = [l l]

    (x, y)
    # nothing
end


""" Plot `LERecordingSolver` via `RecipesBase`."""
@recipe function f(solver::Union{LESolverRecFTLE,
                                 AbstractRenormalizer{<: LESolRecFTLE},
                                 LESolRecFTLE};
                   correlation = true,
                   known_exponents=nothing)
    le_hist, ci_hist = exponents_stat_history(solver)
    dim_lyap = size(le_hist)[1]
    known_exponents = known_exponents == nothing ? [] : known_exponents

    if correlation
        layout := (dim_lyap + 1, 1)
    else
        layout := (dim_lyap, 1)
    end

    ylims = @views [[minimum(le_hist[i, :]),
                     maximum(le_hist[i, :])]
                    for i = 1:dim_lyap]
    for i in 1:min(dim_lyap, length(known_exponents))
        ylims[i][1] = min(ylims[i][1], known_exponents[i])
        ylims[i][2] = max(ylims[i][2], known_exponents[i])
    end
    for i in 1:dim_lyap
        ymin, ymax = ylims[i]
        dy = ymax - ymin
        ylims[i] = [ymin - dy * 0.05,
                    ymax + dy * 0.05]
    end

    for i in 1:dim_lyap
        @series begin
            subplot := i
            label --> ""
            color --> 1
            fillcolor --> 2
            fillalpha --> 0.5
            fillrange := hcat(le_hist[i, 2:end] - ci_hist[i, 2:end],
                              le_hist[i, 2:end] + ci_hist[i, 2:end])
            (2:size(le_hist, 2),
             [le_hist[i, 2:end] le_hist[i, 2:end]])
            # Note: I think it should work w/o "2:end" in principle
            # but it looks like that Plots.jl has some trouble
            # rendering fillrange with NaN.
        end
        if has_convergence_history(solver, i)
            @series begin
                subplot := i
                data = ConvErrorBars(
                    convergence_history(solver),
                    view(le_hist, i, :),
                    i)
                (data,)
            end
        end
        if i <= length(known_exponents)
            @series begin
                subplot := i
                linetype := :hline
                label --> ""
                [known_exponents[i]]
            end
        end
        @series begin
            subplot := i
            label := ""
            ylabel := "LE$i"
            if i == dim_lyap
                xlabel := "Number of orthonormalizations"
            end
            ylim --> ylims[i]
            xlim --> (1, Inf)
            xscale --> :log10
            []
        end
    end

    @series begin
        lags = 0:ceil(Int, sqrt(length(ftle_history(solver))))
        corrs = [autocor(ftle_history(solver, i), lags) for i in 1:dim_lyap]

        subplot := dim_lyap + 1
        ylabel := "Correlation"
        xlabel := "Lag (orthonormalizations)"
        label := ["$i" for i in 1:dim_lyap]
        legend := :top

        (lags, corrs)
    end
end

@recipe function f(demo::LEDemo)
    known_exponents --> demo.example.known_exponents
    demo.solver
end
