using Plots

function list_scripts(dir, ext=".jl")
    paths = []
    for (root, _dirs, files) in walkdir(dir)
        for name in files
            if endswith(name, ext)
                push!(paths, joinpath(root, name))
            end
        end
    end
    return paths
end

for path in list_scripts(joinpath(dirname(@__FILE__), "src/gallery"))
    pngpath = path[1:end-length(".jl")] * ".png"
    if mtime(path) < mtime(pngpath)
        println("Skip running: $path")
        continue
    end
    println("Plotting: $path")
    plt = @time include(path)
    plt = plot(plt, dpi=30)  # plot!(plt, dpi=30) didn't work
    savefig(plt, pngpath)
end
