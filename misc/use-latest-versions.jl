if Base.VERSION < v"0.7-"
    macro info(x)
        :(info($(esc(x))))
    end
else
    using Pkg
end

println("Pkg.free(... all packages ...)")
try
    Pkg.free(collect(Iterators.filter(x -> x != "LyapunovExponents",
                                      keys(Pkg.installed()))))
catch
    # A workaround for Julia 0.7
    for name in keys(Pkg.installed())
        try
            @info "Pkg.free($name)"
            Pkg.free(name)
        end
    end
end
println("Pkg.resolve()")
Pkg.resolve()
println("Pkg.update()")
Pkg.update()
