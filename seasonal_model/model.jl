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
p0 = 0.01
s0 = 0.99                                                           # arbitrary host plant unit
i0 = 0.0                                                     
# encapsulation 
etat0E = @SVector [p0, s0, i0]

# parameters
α = 0.024                                                           # per day
β = 0.04875                                                         # per day per host plant unit
Λ = 0.052                                                           # per day
Θ = 0.04875                                                         # per primary inoculum unit per day
params_elaborate = [α, β, Λ, Θ]


# model for the growing season
function model(u::SVector{3,Float64}, params, t)
    α, β, Λ, Θ = params                                             # unpack the vectors into scalar
    p, s, i = u
    dp = -Λ * p                                                     # dot p
    ds = -Θ * p * s - β * s * i                                     # dot s
    di = Θ * p * s + β * s * i - α * i                              # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

problemE = ODEProblem(model, etat0E, tspan, params_elaborate, saveat=pas_t)
solutionE = solve(problemE)
plot(solutionE, label=["\$P(t)\$" "\$S(t)\$" "\$I(t)\$"], title="Simulation du modèle élaboré")



#=
#etat0C = @SVector [S0, I0]
#α = 0.3698
#β = 0.43
#params_compact = [α, β]
#μ = 0.0072                                                         # per day
μ = 0.005
#π = 1                                                              # arbitrary primary inoculum unit per host plant unit
π = 1.7
λ = 0.2938                                                          # per day
ξ = λ / S0                                                          # per day per host unit
θ = 0.1                                                             # per primary inoculum unit per day
ε = 0.1                                                             # for the elaborate model simulation

function model(u::SVector{2,Float64}, params, t)
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
