#################################################     IMPORTS    ###########################################################################
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector
using Parameters                                                    # for @with_kw
###############################################    STRUCTS AND FUNCTIONS    ################################################################

"""
    struct Growing

Contains the model informations for the growing season,
with few default values.
"""
@with_kw struct Growing
    etat0::SVector{3,Float64}
    params::Vector{Float64}
    tspan::Tuple{Int64,Int64}
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
    all_P::Vector{Vector{Float64}} = []
    all_S::Vector{Vector{Float64}} = []
    all_I::Vector{Vector{Float64}} = []
    all_t::Vector{Vector{Float64}} = []
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
    res.all_P = push!(res.all_P, solutiong[1, :])
    res.all_S = push!(res.all_S, solutiong[2, :])
    res.all_I = push!(res.all_I, solutiong[3, :])
    res.all_t = push!(res.all_t, solutiong.t)

    # useful variables I need below
    s0g = growing.etat0[2]
    π = other.params
    Τ = winter.tspan[2]
    τ = winter.tspan[1]


    # make a first strip
    v1, v2 = [[0, τ]], [[0, τ]]

    # simulation for the rest of the time
    for i in 1:1:years-1

        # complete the vector to make the others stips
        u = [i * Τ, i * Τ + τ]
        push!(v1, u)
        push!(v2, u)

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
        res.all_P = push!(res.all_P, solutionw[1, :])
        res.all_S = push!(res.all_S, solutionw[2, :])
        res.all_I = push!(res.all_I, solutionw[3, :])
        res.all_t = push!(res.all_t, solutionw.t)
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
        res.all_P = push!(res.all_P, solutiong[1, :])
        res.all_S = push!(res.all_S, solutiong[2, :])
        res.all_I = push!(res.all_I, solutiong[3, :])
        res.all_t = push!(res.all_t, solutiong.t)
    end
    
    # convert days into years
    t = res.all_t ./ Τ
    v1 = v1 ./ Τ
    v2 = v2 ./ Τ

    # add stips to plot p1
    p1 = vspan(v1[1],
        color=:lightgray,
        label="growing season")
    p1 = vspan!(v1[2:end],
        color=:lightgray,
        label=false)
    # plot S
    p1 = plot!(t, res.all_S,
        label=false,
        xlims=[0, years],
        ylims=[0, s0g],
        ylabel="\$S(t)\$",
        title="Airborne model",
        legend=:bottomleft,
        c=:black)

    # add stips to plot p2
    p2 = vspan(v2[1],
        color=:lightgray,
        label="growing season")
    p2 = vspan!(v2[2:end],
        color=:lightgray,
        label=false)
    # plot I
    p2 = plot!(t, res.all_I,
        label=false,
        xlims=[0, years],
        ylims=[0, s0g / 3],
        xlabel="Years",
        ylabel="\$I(t)\$",
        legend=:topleft,
        linestyle=:solid,
        c=:black)

    # plot I and P in the same plot, with 2 distincts xaxis
    p2 = plot!(twinx(), t, res.all_P,
        label=false,
        xlims=[0, years],
        ylims=[0, π * s0g / 3],
        ylabel="\$P(t)\$",
        legend=:topright,
        linestyle=:dashdotdot,
        c = :black)

    # subplot S and (P/I)
    plot(p1, p2,
        layout=(2, 1))
    
end    


#######################################################    TEST   ################################################################

t_0 = 0
τ = 184                                                             # growing season length (in days)
Τ = 365                                                             # year duration (in days)
t_transi = τ                                                        # winter season length (in days)
t_fin = Τ
tspang = (t_0, t_transi)
tspanw = (t_transi, t_fin)
temps_simule = 5

# initial conditions
etat0 = @SVector [0.01, 1.0, 0.0]

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
Λ = 0.052                                                           # within-season primary inoculum loss rate per day
Θ = 0.04875                                                         # primary infection rate per primary inoculum unit per day
paramsg = [α, β, Λ, Θ]
μ = 0.0072                                                          # per day

# others parameters
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

growing = Growing(etat0=etat0, params=paramsg, tspan=tspang)
winter = Winter(params=μ, tspan=tspanw)
other = OtherParameters(π)

simule(temps_simule, growing, winter, other)