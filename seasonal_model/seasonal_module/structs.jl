using Parameters, StaticArrays

"""
    `TimeParam` 

Stock all time variables.
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
    `Param` 

This is an abstract type.

Subtypes of `Param` stock bioparameters, model type (compact/elaborate) et the number of states of the model.
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
    `StateParam0` 

This is an abstract type.

Subtypes of `StateParam0` stock model's states at the beginning of a season.
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