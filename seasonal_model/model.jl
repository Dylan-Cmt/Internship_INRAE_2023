# imports 
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 

# time
Τ = 365                                                             # days
τ = 184                                                             # days
#τ = 120

#tspan
t_0 = 0
t_fin = Τ - τ
pas_t = 1
tspan = (t_0, t_fin)

# initial conditions
S0 = 1                                                              # arbitrary host plant unit
I0 = 0
#P0 = 0.01                                                          # for the elaborate model only
# encapsulation 
etat0C = @SVector [S0, I0]
# etat0E = @SVector [P0, S0, I0]
# parameters
#α = 0.024                                                          # per day
α = 0.3698
#β = 0.04875                                                        # per day per host plant unit
β = 0.43
params_compact = [α, β]
Λ = 0.052                                                           # per day
Θ = 0.04875                                                         # per primary inoculum unit per day
# params_elaborate = [α, β, Λ, Θ]
#μ = 0.0072                                                         # per day
μ = 0.005
#π = 1                                                              # arbitrary primary inoculum unit per host plant unit
π = 1.7
λ = 0.2938                                                          # per day
ξ = λ / S0                                                          # per day per host unit
θ = 0.1                                                             # per primary inoculum unit per day
ε = 0.1                                                             # for the elaborate model simulation




# model for the growing season
function model(u::SVector{3,Int64}, params, t)
    α, β, Λ, Θ = params                                             # unpack the vectors into scalar
    P = u[1]
    S = u[2]
    I = u[3]
    dP = -Λ * P                                                     # dot P
    dS = -Θ * P * S - β * S * I                                     # dot S
    dI = Θ * P * S + β * S * I - α * I                              # dot I
    @SVector [dP, dS, dI]                                           # return a new vector
end

problemE = ODEProblem(model, etat0C, tspan, params_compact, saveat=pas_t)
solutionE = solve(problemE)
plot(solutionE, label=["\$S(t)\$" "\$I(t)\$"], title="Simulation du modèle élaboré")

#=
function model(u::SVector{2,Int64}, params, t)
    α, β = params                                                   # unpack the vectors into scalar
    S  = u[1]
    I  = u[2]
    dS = - β * S * I                                                # dot S
    dI = β * S * I - α * I                                          # dot I
    @SVector [dS, dI]                                               # return a new vector
end

problemC = ODEProblem(model, etat0C, tspan, params_compact, saveat=pas_t)
solutionC = solve(problemC)
plot(solutionC, label=["\$S(t)\$" "\$I(t)\$"], title="Simulation du modèle compacte")
=#
