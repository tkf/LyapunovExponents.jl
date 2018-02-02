using RecipesBase

""" Plot `LERecordingSolver` via `RecipesBase`."""
@recipe function f(solver::LERecordingSolver;
                   known_exponents=nothing)
    le_hist = exponents_history(solver)
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
            @view le_hist[i, :]
        end
        if i <= length(known_exponents)
            @series begin
                subplot := i
                linetype := :hline
                label --> ""

                # repeating ylabel/ylim; otherwise they are ignored
                ylabel := "LE$i"
                ylim --> ylims[i]

                [known_exponents[i]]
            end
        end
    end
end

@recipe function f(demo::LEDemo)
    known_exponents --> demo.example.known_exponents
    demo.solver
end
