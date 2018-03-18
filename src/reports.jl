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
end
