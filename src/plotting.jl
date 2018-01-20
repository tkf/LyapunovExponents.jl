using RecipesBase

""" Plot `LERecordingSolver` via `RecipesBase`."""
@recipe function f(solver::LERecordingSolver;
                   known_exponents=nothing)
    exponents_history = solver.exponents_history
    dim_lyap = size(exponents_history)[1]
    layout --> (dim_lyap, 1)
    xscale --> :log10

    ylims = @views [[minimum(exponents_history[i, :]),
                     maximum(exponents_history[i, :])]
                    for i = 1:dim_lyap]
    if known_exponents != nothing
        for i = 1:dim_lyap
            ylims[i][1] = min(ylims[i][1], known_exponents[i])
            ylims[i][2] = max(ylims[i][2], known_exponents[i])
        end
    end
    for i = 1:dim_lyap
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
            @view exponents_history[i, :]
        end
        if known_exponents != nothing
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
