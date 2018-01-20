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
