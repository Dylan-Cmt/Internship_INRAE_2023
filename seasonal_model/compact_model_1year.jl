# imports 
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 
using Parameters                                                    # for @with_kw


"""
    struct Growing

Contains the model informations for the growing season,
with few default values.
"""
@with_kw struct Growing
    etat0::SVector{2,Float64} = @SVector [1.0, 0.0]
    params::Vector{Float64}
    others_params::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    year = 365
    pas = 1
    model::Function = modelg
end


"""
    mutable struct Result

Contains the accumulation of results.
"""
@with_kw mutable struct Result
    all_P::Vector{Vector{Float64}} = []
    all_S::Vector{Vector{Float64}} = []
    all_I::Vector{Vector{Float64}} = []
    all_t::Vector{Vector{Float64}} = []
end


"""
    modelg(u::SVector{3,Float64}, params, t)

It is the model for the growing season.
"""
function modelg(u::SVector{2,Float64}, params, t)
    α, β = params                                                   # unpack the vectors into scalar
    s, i = u
    ds = -β * s * i                                                 # dot s
    di = +β * s * i - α * i                                         # dot i
    @SVector [ds, di]                                               # return a new vector
end

"""
    simule(years, growing::Growing, winter::Winter, res::Result; kwarg...)

Simulates x years of growing seasons.
"""
function simule(years, growing::Growing; kwarg...)

    # Creat a Result type to collect results
    res = Result()
    # collect tspan
    tspan = growing.tspan
    # collect others parameters
    θ, π, μ, λ = growing.others_params

    # inittial conditions
    etat0 = growing.etat0
    @unpack s0, i0 = etat0                                          # susceptible and infected host plant density
    problem = ODEProblem(growing.model, growing.etat0, growing.tspan, growing.params, saveat=growing.pas)
    solution = solve(problem)
    # collect the results
    res.all_S = push!(res.all_S, solutiong[1, :])
    res.all_I = push!(res.all_I, solutiong[2, :])
    res.all_t = push!(res.all_t, solutiong.t)

    # simulation for the rest of the time
    for _ in 1:years-1

        # update tspan of growing season 
        tspan = tspan .+ 365
        # collect the last values to get new initial conditions
        s_fin_g, i_fin_g = last(solution)
        # new initial conditions
        s0g = s0 * exp(-θ * π * exp(-μ(Τ - τ)) * i_fin_g / λ)
        i0g = s0 * (1 - exp(-θ * π * exp(-μ(Τ - τ)) * i_fin_g / λ))
        # encapsulation
        etat0 = @SVector [s0g, i0g]
        # solve the ODE problem for growing season
        problemg = ODEProblem(growing.model, etat0, tspan, growing.params, saveat=growing.pas)
        solutiong = solve(problemg)
        # collect the results
        res.all_S = push!(res.all_S, solutiong[1, :])
        res.all_I = push!(res.all_I, solutiong[2, :])
        res.all_t = push!(res.all_t, solutiong.t)

    end

    return res
end



#=
# plot S
p1 = plot(solutiong.t, [v[1] for v in solutiong.u], label=false,
    xlims  = [0, Τ],
    ylims  = [0, s0_growing + 0.2],
    xlabel = "Year",
    ylabel = "\$S\$")
  
# plot I
p2 = plot(solutiong.t, [v[2] for v in solutiong.u], label=false,
    xlims  = [0, Τ],
    ylims  = [0, s0_growing / 3],
    xlabel = "Year",
    ylabel = "\$I\$")
# plot S et I dans une même fenêtre
plot(p1, p2,
    layout=(2, 1))
title!("Simulation du modèle airborne élaboré", subplot=1)
=#

##############################################    TEST    ################################################################

# time
t_0 = 0
τ = 184
Τ = 365
t_transi = τ
t_fin = Τ

#tspan
pas = 1
tspan = (t_0, t_transi)

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
params = [α, β]
# others parameters
θ = 
π = 
μ = 
λ = 
others_params = [θ, π, μ, λ]

growing = Growing(params=params,others_params = nothing , tspan=(t_0, t_transi), year = Τ)