using Plots

script_dir = "src/gallery/examples"
for name in readdir(script_dir)
    if endswith(name, ".jl")
        path = joinpath(script_dir, name)
        println("Plotting: $path")
        plt = @time include(path)
        savefig(plt, path[1:end-length(".jl")] * ".png")
    end
end
