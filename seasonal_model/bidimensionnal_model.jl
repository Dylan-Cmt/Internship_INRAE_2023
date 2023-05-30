# imports 
using DifferentialEquations                                          # for ODEProblem and solve
using Plots                                                          # for plot
using StaticArrays                                                   # for @SVector 

# time
Τ = 365                                                             # days
τ = 184                                                             # days
#τ = 120

# initial conditions
S0 = 1                                                              # arbitrary host plant unit
I0 = 0
#P0 = 0.01                                                          # for the elaborate model only
# encapsulation 
etat0 = @SVector [S0, I0]

# parameters
#α = 0.024                                                           # per day
α = 0.3698
#β = 0.04875                                                         # per day per host plant unit
β = 0.43
params_compacte = [α, β]
Λ = 0.052                                                           # per day
Θ = 0.04875                                                         # per primary inoculum unit per day
#μ = 0.0072                                                          # per day
μ = 0.005
#π = 1                                                               # arbitrary primary inoculum unit per host plant unit
π = 1.7
λ = 0.2938                                                          # per day
ξ = λ / S0                                                          # per day per host unit
θ = 0.1                                                             # per primary inoculum unit per day
ε = 0.1                                                             # for the elaborate model simulation


#tspan
t_0 = 0.0
t_fin = Τ - τ
pas_t = 0.01
tspan = (t_0, t_fin)

"""
    airborne_sl_compacte_reduced_linearized(u, params, t)

Reduiced and linearized airborne model using the slow-fast argument.
This model is independent of P.
"""
function airborne_sl_compacte_reduced_linearized(u, params, t)
    α, β = params                                                   # unpack the vectors into scalar
    x = u[1]
    y = u[2]
    dx = - β * x * y                                                # dot x
    dy = β * x * y - α * y                                          # dot y
    @SVector [dx, dy]                                               # return a new vector
end


problem = ODEProblem(airborne_sl_compacte_reduced_linearized, etat0, tspan, params_compacte, saveat=pas_t)
solution = solve(problem)
plot(solution,label=["\$S(t)\$" "\$I(t)\$"], title="Simulation du modèle compacte")