# imports 
using DifferentialEquations                                          # for ODEProblem and solve
using Plots                                                          # for plot
using StaticArrays                                                   # for @SVector 
using Parameters

@with_kw struct Modeling
    etat0::Union{SVector{1,Float64},SVector{2,Float64}}              # etat0 is a vector which size is 1 or 2
    p::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas::Float64
    f::Function
end

function simule(mod::Modeling, affiche=false; kwargs...)
    mod_prob = ODEProblem(mod.f, mod.etat0, mod.tspan, mod.p, saveat=mod.pas)
    mod_sol = solve(mod_prob)
    if affiche
        return plot(mod_sol; kwargs...)
    else
        return mod_sol
    end
end

function simule!(mod::Modeling, new_param; kwargs...)
    mod.p[1] = new_param
    mod_prob = ODEProblem(mod.f, mod.etat0, mod.tspan, mod.p, saveat=mod.pas)
    mod_sol = solve(mod_prob)
    plot!(mod_sol; kwargs...)
end

function logistic(u, p, t)
    r, K = p                                        # unpack the vectors into scalar
    x = u[1]
    dx = r * x * (1 - x / K)                        # dot x
    @SVector [dx]                                   # return a new vector
end

etat01 = @SVector [0.1]
r = 1
mod_logistic = Modeling(etat01, [r, 10], (0.0, 10.0), 0.01, logistic)
simule(mod_logistic,true, label="r = $r")

for i in 1.5:0.5:10
    simule!(mod_logistic, i, label="r = $i")
end