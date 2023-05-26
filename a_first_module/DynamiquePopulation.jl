module DynamiquePopulation

# exports
export plot

# imports 
using DifferentialEquations                                          # for ODEProblem and solve
using Plots                                                          # for plot
using StaticArrays                                                   # for @SVector 


struct Modeling
    etat0::Union{SVector{1,Float64},SVector{2,Float64}}             # etat0 is a vector which size is 1 or 2
    p::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas::Float64
    f::Function
end

function ode_solver(mod::Modeling)
    mod_prob = ODEProblem(mod.f, mod.etat0, mod.tspan, mod.p, saveat=mod.pas)
    mod_sol = solve(mod_prob)
    return mod_sol
end

function Plots.plot(mod::Modeling; kwargs...)                       # using ;kwargs... in order to use
    plot(ode_solver(mod); kwargs...)                                # the same arguments as initiale 
end                                                                 # plot fonction

end # end of the module