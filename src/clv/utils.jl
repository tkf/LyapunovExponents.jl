function allocate_array_of_arrays(num, inner_size, inner_type,
                                  constructor = inner_type)
    arrays = Array{inner_type, 1}(num)
    for i in 1:num
        arrays[i] = constructor(inner_size)
    end
    return arrays
end

const UTM = UpperTriangular{Float64, Matrix{Float64}}
make_UTM(args...) = UpperTriangular(Matrix{Float64}(args...))
# TODO: UpperTriangular wastes almost half of memory; fix it
