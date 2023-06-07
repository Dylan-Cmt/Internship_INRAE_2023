#################################################     IMPORTS    ###########################################################################
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector
using Parameters                                                    # for @with_kw
##############################################    PROBLEM INITIALISATION    ################################################################

@with_kw struct Time
    t_0 = 0
    τ                                                               # growing season length (in days)
    Τ = 365                                                         # year duration (in days)
    t_transi = Τ - τ                                                # winter season length (in days)
    t_fin = Τ
end

@with_kw mutable struct Growing
    etat0::SVector{3,Float64} = @SVector [0.0, 0.0, 0.0]
    params::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas::Float64
    model::Function
    solutiong = []
end

@with_kw mutable struct Winter
    etat0::SVector{3,Float64}
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


#=
##############################################    GROWING SEASON: year 1    ################################################################

#tspan
tspang = (t_0, t_transi)

# collect the initial conditions
p_fin_w, s_fin_w, i_fin_w = last(solutionw)

# initial conditions
p0g = p_fin_w                                                       # primary inoculum density
s0g = s_fin_w                                                       # susceptible host plant density
i0g = i_fin_w                                                       # infected host plant density
# encapsulation 
etat0g = @SVector [p0g, s0g, i0g]


problemg = ODEProblem(modelg, etat0g, tspang, params, saveat=pas_t)
solutiong = solve(problemg)

# put solution's values somewhere in order to it plot later
all_P = vcat(all_P, solutiong[1, :])                                # replace last elt by missing for discontinuity
all_S = vcat(all_S, solutiong[2, :])
all_I = vcat(all_I, solutiong[3, :])
all_t = vcat(all_t, solutiong.t)


###########################################    WINTER SEASON: year 1    ################################################################

#tspan
tspanw = (t_transi, t_fin)

# collect growing season data
p_fin_g, s_fin_g, i_fin_g = last(solutiong)

# new initial conditions
p0w = p_fin_g + π * i_fin_g
s0w = 0.0                                                           # arbitrary host plant unit
i0w = 0.0
# encapsulation 
etat0w = @SVector [p0w, s0w, i0w]


problemw = ODEProblem(modelw, etat0w, tspanw, μ, saveat=pas_t)
solutionw = solve(problemw)

# put solution's values somewhere to plot later
all_P = vcat(all_P, missing, solutionw[1, 2:end])
all_S = vcat(all_S, missing, solutionw[2, 2:end-1], missing)        # replace first elt by missing for discontinuity
all_I = vcat(all_I, missing, solutionw[3, 2:end])                   # replace first elt by missing for discontinuity
all_t = vcat(all_t, solutionw.t)


##############################################    GROWING SEASON: year 2    ################################################################

#tspan
tspang = tspang .+ Τ

# collect winter season data
p_fin_w, s_fin_w, i_fin_w = last(solutionw)

# new initial conditions
p0g = p_fin_w
s0g = s0g
i0g = i0g
# encapsulation 
etat0g = @SVector [p0g, s0g, i0g]


problemg = ODEProblem(modelg, etat0g, tspang, params, saveat=pas_t)
solutiong = solve(problemg)

# put solution's values somewhere in order to it plot later
all_P = vcat(all_P, solutiong[1, 1:end-1], missing)                  # replace last elt by missing for discontinuity
all_S = vcat(all_S, solutiong[2, :])
all_I = vcat(all_I, solutiong[3, :])
all_t = vcat(all_t, solutiong.t)
=#

##############################################    TEST   ################################################################

τ = 184                                                             # growing season length (in days)


# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
Λ = 0.052                                                           # within-season primary inoculum loss rate per day
Θ = 0.04875                                                         # primary infection rate per primary inoculum unit per day
params = [α, β, Λ, Θ]
π = 1                                                               # arbitrary primary inoculum unit per host plant unit
μ = 0.0072                                                          # per day

time = Time(τ=184)

winter  = Winter()
growing = Growing()