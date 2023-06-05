# imports 
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 
using Parameters                                                    # for @with_kw

# time
t_0 = 0
τ = 184                                                             # days
Τ = 365                                                             # days
t_transi = Τ - τ
t_fin = Τ

#tspan
pas_t = 1
tspang = (t_0, t_transi)
tspanw = (t_transi, t_fin)


##############################################    GROWING SEASON: year 1    ################################################################


# initial conditions
p0 = 0.01                                                           # primary inoculum density
s0 = 1.0                                                            # susceptible host plant density
i0 = 0.0                                                            # infected host plant density
# encapsulation 
etat0 = @SVector [p0, s0, i0]

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
Λ = 0.052                                                           # within-season primary inoculum loss rate per day
Θ = 0.04875                                                         # primary infection rate per primary inoculum unit per day
π = 1                                                               # arbitrary primary inoculum unit per host plant unit
μ = 0.0072                                                          # per day
params = [α, β, Λ, Θ]

@with_kw struct Mod
    etat0::SVector{3,Float64}
    p::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    tspan::Tuple{Float64,Float64}
    pas::Float64
    f::Function                                                     # modelg or modelw
    season::Bool                                                    # true or false for growing or winter respectively
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

function setInitialCond(mod::Mod)
    if mod.season
        nothing                                                     # initial conditions for growing season
    end
    nothing                                                         # initial conditions for winter season
end

function DifferentialEquations.solve(mod::Mod)
    if mod.season
        mod_prob = ODEProblem(mod.f, mod.etat0, mod.tspan, mod.p, saveat=mod.pas)
        mod_sol = solve(mod_prob)