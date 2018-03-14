packages = [
    let (name, version) = split(line)
        (name, VersionNumber(version))
    end
    for line in split(strip(readstring("REQUIRE")), '\n')[2:end]
]

for (name, version) in packages
    info("Pkg.add($name)")
    Pkg.add(name)  # to make the following work
end
names = map(first, packages)
info("Pkg.free & .update: $(join(names, ", "))")
Pkg.free(names)            # required for Pkg.update
Pkg.update(names...)       # to update METADATA and package repository

for (name, version) in packages
    info("Pkg.pin($name, $version)")
    Pkg.pin(name, version)
end

println("Pkg.resolve()")
Pkg.resolve()
