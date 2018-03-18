report(x) = report(STDOUT, x)

function report(io::IO, sol::LESolution)
    LEs = lyapunov_exponents(sol)

    print_with_color(:blue, io, "Lyapunov Exponents Solution")
    if sol.converged
        print_with_color(:green, io, " (converged)")
    else
        print_with_color(:red, io, " (NOT converged)")
    end
    println(io)

    table = [
        ("#Orth.", sol.num_orth),
        ("#LEs", length(LEs)),
        ("LEs", LEs),
    ]
    for (name, value) in table
        print_with_color(:yellow, io, name)
        print(io, ": ")
        if value isa String
            print(io, value)
        else
            show(IOContext(io, :limit => true), value)
        end
        println(io)
    end

    report(io, sol.convergence)
end

function report(io::IO, convergence::ConvergenceHistory;
                dim_lyap = min(10, length(convergence.errors)))

    print_with_color(:blue, io, "Convergence")
    print(io, " #Orth.=$(convergence.orth[end])")
    if convergence.kinds[end] == UnstableConvError
        print(io, " [unstable]")
    else
        print(io, " [stable]")
    end
    println(io)

    print(io, " "^(length("LE$dim_lyap")), "  ",
          "      error", "   ",
          "  threshold")
    if convergence.kinds[end] == UnstableConvError
        print(io,
              "   variance",
              "   tail cov",
              "  small?")
    end
    println(io)

    for i in 1:dim_lyap
        err = convergence.errors[i][end]
        th = convergence.thresholds[i][end]

        print(io, "LE$i")
        print(io, ": ")
        @printf(io, " %10.5g", err)
        if err < th
            print_with_color(:green, io, " < ")
        else
            print_with_color(:red, io, " > ")
        end
        @printf(io, " %10.5g", th)

        if convergence.kinds[end] == UnstableConvError
            detail = convergence.details[i][end]
            @printf(io, " %10.5g", detail.var)
            @printf(io, " %10.5g", detail.tail_cov)
            print(io, "  ")
            if detail.tail_ok
                print_with_color(:green, io, "yes")
            else
                print_with_color(:red, io, "no")
            end
        end

        println(io)
    end
end
