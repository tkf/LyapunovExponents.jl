module Examples

include("example_base.jl")

include("arnold_cat_map.jl")
include("bakers_map.jl")
include("henon_map.jl")
include("linz_sprott_99.jl")
include("lorenz_63.jl")
include("lozi_map.jl")
include("rnn.jl")
include("standard_map.jl")
include("van_der_pol.jl")

end

using .Examples.ExampleBase: LEDemo, LEExample

# Make LyapunovExponents.lorenz_63 etc. available:
using .Examples.ArnoldCatMap
using .Examples.BakersMap
using .Examples.HenonMap
using .Examples.LinzSprott99
using .Examples.Lorenz63
using .Examples.LoziMap
using .Examples.RNN
using .Examples.StandardMap
using .Examples.VanDerPol

const DEMOS = [
    lorenz_63,
    linz_sprott_99,
    van_der_pol,
    beer_95,
    henon_map,
    standard_map,
    bakers_map,
    arnold_cat_map,
    lozi_map,
]
