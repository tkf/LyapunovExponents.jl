if Base.VERSION < v"0.7-"
    _readstring = readstring
else
    _readstring = f -> read(f, String)
    using Pkg
end

packages = Pkg.installed()
for line in split(strip(_readstring("REQUIRE")), '\n')[2:end]
    name = split(line)[1]
    println(name, "\t", get(packages, name, ""))
end
