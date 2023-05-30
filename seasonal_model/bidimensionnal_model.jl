# imports 
using DifferentialEquations                                          # for ODEProblem and solve
using Plots                                                          # for plot
using StaticArrays                                                   # for @SVector 

# time
Τ = 365                                                             # days
τ = 184                                                             # days


# initial conditions
S0 = 1 # arbitrary host plant unit
I0 = 0
#P0 = 0
# encapsulation 
etat0 = @SVector [S0, I0]

# parameters
α = 0.024                                                           # per day
β = 0.04875                                                         # per day per host plant unit
params_compacte = [α, β]
Λ = 0.052                                                           # per day
Θ = 0.04875                                                         # per primary inoculum unit per day
ε = 0.1
λ = ε * Λ
θ = ε * Θ
π = 1                                                               # arbitrary primary inoculum unit per host plant unit
μ = 0.0072                                                          # per day

ξ = λ / S0                                                          # per day per host unit


#tspan
t_0 = 0.0
t_fin = Τ - τ
pas_t = 0.01
tspan = (t_0, t_fin)


function model_compacte(u, params, t)
    α, β = params                                                   # unpack the vectors into scalar
    x = u[1]
    y = u[2]
    dx = - β * x * y                                                # dot x
    dy = β * x * y - α * y                                          # dot y
    @SVector [dx, dy]                                               # return a new vector
end


problem = ODEProblem(model_compacte, etat0, tspan, params_compacte, saveat=pas_t)
solution = solve(problem)
plot(solution)