using Parameters, StaticArrays

"""
`TimeParam` stock toutes les variables relatives au temps
"""
@with_kw struct TimeParam
    T::Float64 = 365
    @assert T > 0
    @assert T <= 365 # year length
    τ::Float64 = 184
    @assert τ <= T # crop season length
    Δt::Float64 = 0.1
    @assert Δt > 0 # step
    tspang::Tuple{Float64,Float64} = (0, τ) # tspan for growing
    tspanw::Tuple{Float64,Float64} = (τ, T) # tspan for winter	
end

"""
`Param` stock les paramètres biologiques, le type de modèle (compact/élaboré) et le nombre d'états du modèle.

Les différentes struct ont été trié à l'aide d'abstract types, de sorte à pouvoir utiliser des types plus généraux dans nos fonctions. 

Par exemple, `ParamAirborneElaborate1Strain` et `ParamSoilborneElaborate1Strain` ont été regroupé dans `Elaborate1Strain` car ils partagent la méthode `WinterSeason`.
"""
abstract type Param end


abstract type Elaborate1Strain <: Param end
abstract type Elaborate2Strains <: Param end
abstract type Compact1Strain <: Param end
abstract type Compact2Strains <: Param end


@with_kw struct ParamSoilborneCompact1Strain{T<:Float64} <: Compact1Strain
    # growing season parameters
    α::T = 0.024
    @assert α > 0 # Infected host plants removal rate
    β::T = 0.04875
    @assert β > 0 # Secondary infection rate

    # convertion into primary inoculum parameter
    θ::T = 0.04875
    @assert θ > 0 # Primary infection rate
    Π::T = 1.0
    @assert Π > 0 # Conversion rate from I to P (at the end of the season)
    μ::T = 0.0072
    @assert μ > 0
    λ::T = 0.0072
    @assert λ > 0
    n::T = 1.0
    @assert n >= 0 # Initial plant density

    # type of model and state number 
    isElaborate = false
    statelength = 2
end

@with_kw struct ParamAirborneElaborate1Strain{T<:Float64} <: Elaborate1Strain
    # growing season parameters
    Λ::T = 0.052
    @assert Λ > 0 # Primary inoculum density independent depletion rate
    Θ::T = 0.04875
    @assert Θ > 0 # Primary infection rate
    α::T = 0.024
    @assert α > 0 # Infected host plants removal rate
    β::T = 0.04875
    @assert β > 0 # Secondary infection rate

    # convertion into primary inoculum parameter
    Π::T = 1.0
    @assert Π > 0 # Conversion rate from I to P (at the end of the season)

    # winter-specific mortality parameter
    μ::T = 0.0072
    @assert μ > 0 # Winter season mortality rate of primary inoculum

    # new susceptible host plant density parameter 
    n::T = 1.0
    @assert n >= 0 # Initial plant density

    # type of model and state number 
    isElaborate = true
    statelength = 3
end

"""
`StateParam0` stock les états du modèle au début d'une saison.
"""
abstract type StateParam0 end

@with_kw struct StateCompact <: StateParam0
    S0::Float64 = 0.99
    @assert S0 >= 0
    I0::Float64 = 0.01
    @assert I0 >= 0
    @assert S0 + I0 <= 1
    State0 = @SVector [S0, I0]
end

@with_kw struct StateElaborate <: StateParam0
    P0::Float64 = 0.01
    @assert P0 >= 0
    S0::Float64 = 0.99
    @assert S0 >= 0
    I0::Float64 = 0.00
    @assert I0 >= 0
    @assert S0 + I0 <= 1
    State0 = @SVector [P0, S0, I0]
end

