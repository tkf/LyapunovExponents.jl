if Base.VERSION >= v"0.7"
    using Pkg
end

println("Pkg.free(... all packages ...)")
Pkg.free(collect(Iterators.filter(x -> x != "LyapunovExponents",
                                  keys(Pkg.installed()))))
println("Pkg.resolve()")
Pkg.resolve()
println("Pkg.update()")
Pkg.update()
