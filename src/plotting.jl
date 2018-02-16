using RecipesBase

""" Plot `LERecordingSolver` via `RecipesBase`."""
@recipe function f(solver::Union{LESolverRecFTLE,
                                 AbstractRenormalizer{<: LESolRecFTLE},
                                 LESolRecFTLE};
                   known_exponents=nothing)
    le_hist, ci_hist = exponents_stat_history(solver)
    dim_lyap = size(le_hist)[1]
    known_exponents = known_exponents == nothing ? [] : known_exponents
    layout --> (dim_lyap, 1)
    xscale --> :log10

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
            ylabel := "LE$i"
            ylim --> ylims[i]
            xlim --> (1, Inf)
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
        if i <= length(known_exponents)
            @series begin
                subplot := i
                linetype := :hline
                label --> ""

                # repeating ylabel/ylim; otherwise they are ignored
                ylabel := "LE$i"
                ylim --> ylims[i]
                xlim --> (1, Inf)

                [known_exponents[i]]
            end
        end
    end
end

@recipe function f(demo::LEDemo)
    known_exponents --> demo.example.known_exponents
    demo.solver
end
