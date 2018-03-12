packages = [
    let (name, version) = split(line)
        (name, VersionNumber(version))
    end
    for line in split(strip(readstring("REQUIRE")), '\n')[2:end]
]

for (name, version) in packages
    info("Pkg.add($name)")
    Pkg.add(name)  # to make the following work
    info("Pkg.free($name)")
    Pkg.free(name)  # required for Pkg.update
    info("Pkg.update($name)")
    Pkg.update(name)  # to update METADATA and package repository
    info("Pkg.pin($name, $version)")
    Pkg.pin(name, version)
end
println("Pkg.resolve()")
Pkg.resolve()
