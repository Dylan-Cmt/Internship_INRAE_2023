#################################################     IMPORTS    ###########################################################################
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector
using Parameters                                                    # for @with_kw
##############################################    PROBLEM INITIALISATION    ################################################################


@with_kw mutable struct Growing
    etat0::SVector{3,Float64} = @SVector [0.0, 0.0, 0.0]
    params::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas::Float64
    model::Function
    solutiong = []
end

@with_kw mutable struct Winter
    etat0::SVector{3,Float64} = @SVector [0.0, 0.0, 0.0]
    μ::Float64
    tspan::Tuple{Float64,Float64}
    pas::Float64
    model::Function
    solutionw = [[0.01, 1.0, 0.0]]                                  # contains initial state of the problem
end

# accumulation of results (also type missing to plot a discontinuity)
@with_kw mutable struct Result
    all_P::Vector{Union{Missing,Float64}} = []
    all_S::Vector{Union{Missing,Float64}} = []
    all_I::Vector{Union{Missing,Float64}} = []
    all_t::Vector{Union{Missing,Float64}} = []
end

# model for the growing season
function modelg(u::SVector{3,Float64}, params, t)
    α, β, Λ, Θ = params                                             # unpack the vectors into scalar
    p, s, i = u
    dp = -Λ * p                                                     # dot p
    ds = -Θ * p * s - β * s * i                                     # dot s
    di = Θ * p * s + β * s * i - α * i                              # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

# model for the winter season
function modelw(u::SVector{3,Float64}, params, t)
    μ = params                                                      # unpack the vectors into scalar
    p, s, i = u
    dp = -μ * p                                                     # dot p
    ds = 0                                                          # dot s
    di = 0                                                          # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

function simule(years, growing::Growing, winter::Winter, res::Result, plot=true; kwarg...)
    # simulate n years
    for i in 1:years
        # growing season
        growing.tspan = growing.tspan .+ (i-1)*365
        p_fin_w, s_fin_w, i_fin_w = last(winter.solutionw)
        p0g = p_fin_w
        s0g = s_fin_w
        i0g = i_fin_w
        growing.etat0 = @SVector [p0g, s0g, i0g]
        problemg = ODEProblem(growing.modelg, growing.etat0, growing.tspan, growing.params, saveat=growing.pas_t)
        growing.solutiong = solve(problemg)

        res.all_P = vcat(res.all_P, growing.solutiong[1, :])
        res.all_S = vcat(res.all_S, growing.solutiong[2, :])
        res.all_I = vcat(res.all_I, growing.solutiong[3, :])
        res.all_t = vcat(res.all_t, growing.solutiong.t)


        # winter season
        winter.tspan = winter.tspan .+ (i - 1) * 365
        p_fin_g, s_fin_g, i_fin_g = last(winter.solutionw)
        p0w = p_fin_g
        s0w = s_fin_g
        i0w = i_fin_g
        winter.etat0 = @SVector [p0w, s0w, i0w]
        problemw = ODEProblem(winter.modelw, winter.etat0, winter.tspan, winter.μ, saveat=winter.pas_t)
        winter.solutionw = solve(problemw)

        res.all_P = vcat(res.all_P, missing, winter.solutionw[1, 2:end])
        res.all_S = vcat(res.all_S, missing, winter.solutionw[2, 2:end-1], missing)
        res.all_I = vcat(res.all_I, missing, winter.solutionw[3, 2:end])
        res.all_t = vcat(res.all_t, solutionw.t)
    end
    if plot
    # plot the result
        # oublies pas de faire Plots.plot(...)
    end
end    


##############################################    TEST   ################################################################

t_0 = 0
τ = 184                                                             # growing season length (in days)
Τ = 365                                                             # year duration (in days)
t_transi = Τ - τ                                                    # winter season length (in days)
t_fin = Τ

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
Λ = 0.052                                                           # within-season primary inoculum loss rate per day
Θ = 0.04875                                                         # primary infection rate per primary inoculum unit per day
params = [α, β, Λ, Θ]
π = 1                                                               # arbitrary primary inoculum unit per host plant unit
μ = 0.0072                                                          # per day


growing = Growing(params=params, tspan=(t_0, t_transi), pas=1, model=modelg)
winter = Winter(μ=μ, tspan=(t_transi, t_fin), pas=1, model=modelw)
res = Result()
