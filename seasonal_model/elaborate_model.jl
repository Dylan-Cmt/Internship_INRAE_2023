#################################################     IMPORTS    ###########################################################################
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector
using Parameters                                                    # for @with_kw
##############################################    PROBLEM INITIALISATION    ################################################################

"""
    struct Growing

Contains the model informations for the growing season,
with few default values.
"""
@with_kw struct Growing
    etat0::SVector{3,Float64}
    params::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas = 1
    model::Function = modelg
end

"""
    struct Winter

Contains the model informations for the winter season,
with few default values.
"""
@with_kw struct Winter
    params::Union{Float64,Vector{Float64}}
    tspan::Tuple{Int64,Int64}
    pas = 1
    model::Function = modelw
end

"""
    struct OtherParameters

Contains other parameters.
"""
mutable struct OtherParameters
    params::Union{Int64,Float64,Vector{Float64}}
end

"""
    mutable struct Result

Contains the accumulation of results.
Can also contains Missing type to include discontinuity.
"""
@with_kw mutable struct Result
    all_P::Vector{Union{Missing,Float64}} = []
    all_S::Vector{Union{Missing,Float64}} = []
    all_I::Vector{Union{Missing,Float64}} = []
    all_t::Vector{Union{Missing,Float64}} = []
end

"""
    modelg(u::SVector{3,Float64}, params, t)

It is the model for the growing season.
"""
function modelg(u::SVector{3,Float64}, params, t)
    α, β, Λ, Θ = params                                             # unpack the vectors into scalar
    p, s, i    = u
    dp = -Λ * p                                                     # dot p
    ds = -Θ * p * s - β * s * i                                     # dot s
    di = Θ * p * s + β * s * i - α * i                              # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

"""
    modelw(u::SVector{3,Float64}, params, t)

It is the model for the winter season.
"""
function modelw(u::SVector{3,Float64}, params, t)
    μ       = params                                                # unpack the vectors into scalar
    p, s, i = u
    dp = -μ * p                                                     # dot p
    ds = 0                                                          # dot s
    di = 0                                                          # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

"""
    simule(years, growing::Growing, winter::Winter, res::Result; kwarg...)

Simulates x years of alternation between growing and winter seasons.
"""
function simule(years, growing::Growing, winter::Winter, other::OtherParameters; kwarg...)
    
    # Creat a Result type to collect results
    res = Result()

    # collect the tspans
    tspang = growing.tspan
    tspanw = winter.tspan

    # GROWING SEASON
    # solve the ODE problem for a first growing season
    problemg  = ODEProblem(growing.model, growing.etat0, tspang, growing.params, saveat=growing.pas)
    solutiong = solve(problemg)
    # collect the results
    res.all_P = vcat(res.all_P, solutiong[1, :])
    res.all_S = vcat(res.all_S, solutiong[2, :])
    res.all_I = vcat(res.all_I, solutiong[3, :])
    res.all_t = vcat(res.all_t, solutiong.t)

    # this variable will never change so let's put it before the loop
    s0g = growing.etat0[2]
    π = other.params

    # simulation for the rest of the time
    for _ in 1:years-1

        # WINTER SEASON
        # collect the last values to get new initial conditions
        p_fin_g, s_fin_g, i_fin_g = last(solutiong)
        # new initial conditions
        p0w = p_fin_g + π * i_fin_g
        s0w = 0.0
        i0w = 0.0
        # encapsulation 
        etat0 = @SVector [p0w, s0w, i0w]
        # solve the ODE problem for winter season
        problemw = ODEProblem(winter.model, etat0, tspanw, winter.params, saveat=winter.pas)
        solutionw = solve(problemw)
        # collect the results
        res.all_P = vcat(res.all_P, missing, solutionw[1, 2:end])
        res.all_S = vcat(res.all_S, missing, solutionw[2, 2:end-1], missing)
        res.all_I = vcat(res.all_I, missing, solutionw[3, 2:end])
        res.all_t = vcat(res.all_t, solutionw.t)
        # update tspan of winter season for the next year
        tspanw = tspanw .+ 365

        # GROWING SEASON
        # update tspan of growing season 
        tspang = tspang .+ 365
        # collect the last values to get new initial conditions
        p_fin_w, s_fin_w, i_fin_w = last(solutionw)
        # new initial conditions
        p0g = p_fin_w
        i0g = 0.0
        # encapsulation
        etat0 = @SVector [p0g, s0g, i0g]
        # solve the ODE problem for growing season
        problemg = ODEProblem(growing.model, etat0, tspang, growing.params, saveat=growing.pas)
        solutiong = solve(problemg)
        # collect the results
        res.all_P = vcat(res.all_P, solutiong[1, :])
        res.all_S = vcat(res.all_S, solutiong[2, :])
        res.all_I = vcat(res.all_I, solutiong[3, :])
        res.all_t = vcat(res.all_t, solutiong.t) 
    end
    
    # convert days into years
    t = res.all_t ./ winter.tspan[2]

    # plot I
    p1 = Plots.plot(t, res.all_I,
        label="\$I\$",
        legend=:topleft,
        c=:red,
        xlabel="Years",
        ylabel="\$I(t)\$",
        linestyle=:solid,
        ylims=[0, s0g / 3])

    # plot I and P in the same plot
    p1 = Plots.plot!(twinx(), t, res.all_P,
        c=:black,
        label="\$P\$",
        legend=:topright,
        ylabel="\$P(t)\$",
        linestyle=:dashdotdot,
        ylims=[0, π * s0g / 3])

    # plot S
    p2 = Plots.plot(t, res.all_S,
        ylims=[0, s0g],
        label=false,
        ylabel="\$S(t)\$",
        title="Airborne model")

    # subplot S and (P/I)
    Plots.plot(p2, p1,
        layout=(2, 1),
        xlims=[0, years])
    
end    


#######################################################    TEST   ################################################################

t_0 = 0
τ = 184                                                             # growing season length (in days)
Τ = 365                                                             # year duration (in days)
t_transi = τ                                                        # winter season length (in days)
t_fin = Τ

etat0 = @SVector [0.01, 1.0, 0.0]

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
Λ = 0.052                                                           # within-season primary inoculum loss rate per day
Θ = 0.04875                                                         # primary infection rate per primary inoculum unit per day
paramsg = [α, β, Λ, Θ]
μ = 0.0072                                                          # per day
π = 1                                                               # arbitrary primary inoculum unit per host plant unit
#=
# graph Fig. 3
τ = 120
α = 0.3698                                                          
β = 0.43                                                            
λ = 0.2938
ϵ = 0.1
Λ = λ/ϵ                                                             
θ =0.1
Θ = θ/ϵ                                                             
paramsg = [α, β, Λ, Θ]
μ = 0.005                                                           
π = 1.7                                                             
=#

growing = Growing(etat0, paramsg, (t_0, t_transi))
winter = Winter(params=μ, tspan=(t_transi, t_fin))
other = OtherParameters(π)

simule(5, growing, winter, other)