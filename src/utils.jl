using DifferentialEquations: DEProblem
using ProgressMeter

macro showprogress_if(pred, args...)
    quote
        if $(esc(pred))
            @showprogress $(map(esc, args)...)
        else
            $(map(esc, args)...)
        end
    end
end

objname(f) = rsplit(string(f), "."; limit=2)[end]

function is_semi_unitary(U, error=1e-7)
    n, m = size(U)
    I = n > m ? U' * U : U * U'
    diff = abs.(I .- eye(I))
    sum(diff) / length(diff) <= error
end

function default_Q0(T::DataType, dim_phase, dim_lyap)
    Q0 = zeros(T, dim_phase, dim_lyap)
    for i in 1:dim_lyap
        Q0[1:dim_phase - i + 1, i] = 1
    end
    return qr(Q0)[1]
end

default_Q0(array::AbstractArray, dim_phase, dim_lyap) =
    default_Q0(eltype(array), dim_phase, dim_lyap)

default_Q0(prob::DEProblem, dim_phase, dim_lyap) =
    default_Q0(prob.u0, dim_phase, dim_lyap)

default_Q0(x, dim) = default_Q0(x, dim, dim)
