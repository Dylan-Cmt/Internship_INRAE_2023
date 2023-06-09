#################################################     IMPORTS    ###########################################################################
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector
using Parameters                                                    # for @with_kw

###############################################    STRUCTS AND FUNCTIONS    ################################################################

"""
    struct Growing

Contains the model informations for the growing season.
"""
@with_kw struct Growing
    etat0::SVector{5,Float64} = @SVector [0.01, 0.01, 1.0, 0.0, 0.0]
    params::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas = 1
    modelg::Function
end

"""
    struct Winter

Contains the model informations for the winter season,
with few default values.
"""
@with_kw struct Winter
    params::Union{Float64,Vector{Float64}}
    tspan::Tuple{Float64,Float64}
    pas = 1
    modelw::Function
end

"""
    struct OtherParameters

Contains other parameters.
"""
@with_kw struct OtherParameters
    params::Union{Float64,Vector{Float64}}
end

"""
    mutable struct Result

Contains the accumulation of results.
Can also contains Missing type to include discontinuity.
"""
@with_kw mutable struct Result
    p1::Vector{Union{Missing,Float64}} = []
    p2::Vector{Union{Missing,Float64}} = []
    s::Vector{Union{Missing,Float64}}  = []
    i1::Vector{Union{Missing,Float64}} = []
    i2::Vector{Union{Missing,Float64}} = []
    t::Vector{Union{Missing,Float64}}  = []
end

"""
    modelg(u::SVector{3,Float64}, params, t)

It is the model for the growing season.
"""
function modelg(u::SVector{3,Float64}, params, t)
    α, β1, β2, Λ, Θ   = params                                      # unpack the vectors into scalar
    p1, p2, s, i1, i2 = u
    dp1 = -Λ * p1                                                   # dot p1
    dp2 = -Λ * p2                                                   # dot p2
    ds  = -Θ * p1 * s + β1 * s * i1   -Θ * p2 * s + β2 * s * i2     # dot s
    di1 = Θ * p1 * s + β1 * s * i1 - α * i1                         # dot i1
    di2 = Θ * p2* s + β2 * s * i2 - α * i1  2                       # dot i2
    @SVector [dp1, dp2, ds, di1, di2]                                  # return a new vector
end

"""
    modelw(u::SVector{3,Float64}, params, t)

It is the model for the winter season.
"""
function modelw(u::SVector{3,Float64}, params, t)
    μ1, μ2            = params                                      # unpack the vectors into scalar
    p1, p2, s, i1, i2 = u
    dp1 = -μ1 * p1                                                  # dot p1
    dp2 = -μ2 * p2                                                  # dot p2
    ds  = 0                                                         # dot s
    di1 = 0                                                         # dot i1
    di2 = 0                                                         # dot i2
    @SVector [dp1, dp2, ds, di1, di2]                               # return a new vector
end

"""
    simule(years,model::Function, growing::Growing, winter::Winter)

Simulates x years of alternation between growing and winter seasons.
"""
function simule(years, growing::Growing, winter::Winter)
    # Creat a Result type to collect results
    res = Result()

    # collect the tspans
    tspang = growing.tspan
    tspanw = winter.tspan

    # GROWING SEASON
    # solve the ODE problem for a first growing season
    problemg  = ODEProblem(growing.modelg, growing.etat0, tspang, growing.params, saveat=growing.pas)
    solutiong = solve(problemg)
    # collect the results
    res.p1 = vcat(res.p1, solutiong[1, :])
    res.p1 = vcat(res.p1, solutiong[2, :])
    res.s  = vcat(res.s, solutiong[3, :])
    res.i1 = vcat(res.i1, solutiong[4, :])
    res.i2 = vcat(res.i2, solutiong[5, :])
    res.t  = vcat(res.t, solutiong.t)    
end

#######################################################    TEST   ################################################################

t_0 = 0
τ   = 184                                                           # growing season length (in days)
Τ   = 365                                                           # year duration (in days)
t_transi = τ                                                        # winter season length (in days)
t_fin    = Τ

# parameters
α  = 0.024                                                          # infected host plants removal rate per day
β1 = 0.04875                                                        # secondary infection rate per day per host plant unit
β2 = 0.04875                                                        # secondary infection rate per day per host plant unit
Λ  = 0.052                                                          # within-season primary inoculum loss rate per day
Θ  = 0.04875                                                        # primary infection rate per primary inoculum unit per day
paramsg = [α, β1, β2, Λ, Θ]
μ1 = 0.0072                                                         # per day
μ2 = 0.0072                                                         # per day
paramsw = [ μ1, μ2]
π  = 1                                                              # arbitrary primary inoculum unit per host plant unit

growing1 = Growing(params=paramsg, tspan=(t_0, t_transi))
winter1  = Winter(params=paramsw, tspan=(t_transi, t_fin))
other1   = OtherParameters(π)
