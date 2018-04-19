if Base.VERSION >= v"0.7"
    using Pkg
end

dev_packages = ["DiffEqBase", "OrdinaryDiffEq", "StochasticDiffEq",
                "OnlineStats", "OnlineStatsBase"]
dont_free = vcat(dev_packages, ["LyapunovExponents"])

info("Pkg.free(... all packages ...)")
Pkg.free(collect(Iterators.filter(x -> ! (x in dont_free),
                                  keys(Pkg.installed()))))
for name in dev_packages
    info("Pkg.add($name)")
    Pkg.add(name)  # to make the following work
    info("Pkg.checkout($name)")
    Pkg.checkout(name)
end
info("Pkg.resolve()")
Pkg.resolve()
info("Pkg.update()")
Pkg.update()
