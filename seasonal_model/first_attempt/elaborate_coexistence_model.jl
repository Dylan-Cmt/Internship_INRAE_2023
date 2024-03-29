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
    etat0::SVector{5,Float64}
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
    params::Union{Int64, Float64,Vector{Float64}}
end

"""
    mutable struct Result

Contains the accumulation of results.
Can also contains Missing type to include discontinuity.
"""
@with_kw mutable struct Result
    p1::Vector{Vector{Float64}} = []
    p2::Vector{Vector{Float64}} = []
    s::Vector{Vector{Float64}} = []
    i1::Vector{Vector{Float64}} = []
    i2::Vector{Vector{Float64}} = []
    t::Vector{Vector{Float64}} = []
end

"""
    modelg(u::SVector{3,Float64}, params, t)

It is the model for the growing season.
"""
function modelg(u::SVector{5,Float64}, params, t)
    α, β1, β2, Λ, Θ   = params                                      # unpack the vectors into scalar
    p1, p2, s, i1, i2 = u
    dp1 = -Λ * p1                                                   # dot p1
    dp2 = -Λ * p2                                                   # dot p2
    ds  = -Θ * p1 * s + β1 * s * i1   -Θ * p2 * s + β2 * s * i2     # dot s
    di1 = Θ * p1 * s + β1 * s * i1 - α * i1                         # dot i1
    di2 = Θ * p2* s + β2 * s * i2 - α * i1                          # dot i2
    @SVector [dp1, dp2, ds, di1, di2]                               # return a new vector
end

"""
    modelw(u::SVector{3,Float64}, params, t)

It is the model for the winter season.
"""
function modelw(u::SVector{5,Float64}, params, t)
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
function simule(years, growing::Growing, winter::Winter, other::OtherParameters)

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
    res.p1 = push!(res.p1, solutiong[1, :])
    res.p2 = push!(res.p2, solutiong[2, :])
    res.s = push!(res.s, solutiong[3, :])
    res.i1 = push!(res.i1, solutiong[4, :])
    res.i2 = push!(res.i2, solutiong[5, :])
    res.t = push!(res.t, solutiong.t)
    
    π  = other.params
    s0 = growing.etat0[3]

    # simulation for the rest of the time
    for _ in 1:1:years-1

        # WINTER SEASON
        # collect the last values to get new initial conditions
        p1_fin_g, p2_fin_g, s_fin_g, i1_fin_g, i2_fin_g = last(solutiong)
        # new initial conditions
        p01w = p1_fin_g + π * i1_fin_g
        p02w = p2_fin_g + π * i2_fin_g
        s0w  = 0.0
        i01w = 0.0
        i02w = 0.0
        # encapsulation 
        etat0 = @SVector [p01w, p02w, s0w, i01w, i02w]
        # solve the ODE problem for winter season
        problemw = ODEProblem(winter.model, etat0, tspanw, winter.params, saveat=winter.pas)
        solutionw = solve(problemw)
        # collect the results
        res.p1 = push!(res.p1, solutionw[1, :])
        res.p2 = push!(res.p2, solutionw[2, :])
        res.s = push!(res.s, solutionw[3, :])
        res.i1 = push!(res.i1, solutionw[4, :])
        res.i2 = push!(res.i2, solutionw[5, :])
        res.t = push!(res.t, solutionw.t)
        # update tspan of winter season for the next year
        tspanw = tspanw .+ 365

        # GROWING SEASON
        # update tspan of growing season 
        tspang = tspang .+ 365
        # collect the last values to get new initial conditions
        p1_fin_w, p2_fin_w, s_fin_w, i1_fin_w, i2_fin_w = last(solutionw)
        # new initial conditions
        p01g = p1_fin_g
        p02g = p2_fin_g
        s0g  = s0
        i01g = 0.0
        i02g = 0.0
        # encapsulation
        etat0 = @SVector [p01g, p02g, s0g, i01g, i02g]
        # solve the ODE problem for growing season
        problemg = ODEProblem(growing.model, etat0, tspang, growing.params, saveat=growing.pas)
        solutiong = solve(problemg)
        # collect the results
        res.p1 = push!(res.p1, solutiong[1, :])
        res.p2 = push!(res.p2, solutiong[2, :])
        res.s = push!(res.s, solutiong[3, :])
        res.i1 = push!(res.i1, solutiong[4, :])
        res.i2 = push!(res.i2, solutiong[5, :])
        res.t = push!(res.t, solutiong.t)
    end
    # plot
    t = res.t ./ Τ

    # plot S
    plt1 = plot(t, res.s,
        label=false,
        xlims=[0, years],
        ylims=[0, s0],
        ylabel="\$S(t)\$",
        title="Airborne model",
        legend=:bottomleft,
        c=:black)

    # plot p1
    plt2 = plot(t, res.p1,
        label=false,
        xlims=[0, years],
        ylims=[0, s0 / 3],
        xlabel="Years",
        ylabel="\$P1(t)\$",
        legend=:topleft,
        linestyle=:solid,
        c=:black)
    

    # plot I and P in the same plot, with 2 distincts xaxis
    plt2 = plot!(twinx(), t, res.i1,
        label=false,
        xlims=[0, years],
        ylims=[0, π * s0 / 3],
        ylabel="\$I1(t)\$",
        legend=:topright,
        linestyle=:dashdotdot,
        c=:black)

    # subplot S and (P/I)
    plot(plt1, plt2,
        layout=(2, 1))
end

#######################################################    TEST   ################################################################

# time
t_0 = 0
τ   = 165                                                           # growing season length (in days)
Τ   = 365                                                           # year duration (in days)
t_transi = τ                                                        # winter season length (in days)
t_fin    = Τ
tspang = (t_0, t_transi)
tspanw = (t_transi, t_fin)
temps_simule = 2

# initial conditions
etat0 = @SVector [0.01, 0.01, 1.0, 0.0, 0.0]

# parameters
α  = 0.005                                                          # infected host plants removal rate per day
β1 = 0.01                                                           # secondary infection rate per day per host plant unit
β2 = 0.035                                                          # secondary infection rate per day per host plant unit
Λ  = 0.01                                                           # within-season primary inoculum loss rate per day
Θ  = 0.05                                                           # primary infection rate per primary inoculum unit per day
paramsg = [α, β1, β2, Λ, Θ]
μ1 = 0.0025                                                         # per day
μ2 = 0.0068                                                         # per day
paramsw = [μ1, μ2]

# others parameters
π  = 1                                                              # arbitrary primary inoculum unit per host plant unit

growing = Growing(etat0=etat0, params=paramsg, tspan=tspang)
winter = Winter(params=paramsw, tspan=tspanw)
other   = OtherParameters(π)

resultat = simule(temps_simule, growing, winter, other)
