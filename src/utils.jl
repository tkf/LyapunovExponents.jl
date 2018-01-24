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

function is_semi_unitary(U, error=1e-7)
    n, m = size(U)
    I = n > m ? U' * U : U * U'
    diff = abs.(I .- eye(I))
    sum(diff) / length(diff) <= error
end
