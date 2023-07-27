include("Seasonal.jl")
using .Seasonal

# time parameter
tp = TimeParam()
# elaborate 1 strain model
paramE = ParamAirborneElaborate1Strain()
spE = StateElaborate()
# compact 1 strain model
paramC = ParamSoilborneCompact1Strain()
spC = StateCompact()


affiche(4, spE, paramE)