packages = [
    let (name, version) = split(line)
        (name, VersionNumber(version))
    end
    for line in split(strip(readstring("REQUIRE")), '\n')[2:end]
]

for (name, version) in packages
    println("Pkg.pin($name, $version)")
    Pkg.pin(name, version)
end
println("Pkg.resolve()")
Pkg.resolve()
