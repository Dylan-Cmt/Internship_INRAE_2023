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
    etat0::SVector{2,Float64} = @SVector [1.0, 0.01]
    params::Vector{Float64}
    others_params::Vector{Float64}
    tspan::Tuple{Int64,Int64}
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
    simule(years, growing::Growing)

Simulates x years of growing seasons.
"""
function simule(years, growing::Growing)

    # Creat a Result type to collect results
    res = Result()

    # collect tspan
    tspan = growing.tspan

    # collect others parameters
    θ, π, μ, λ = growing.others_params

    # initial condition
    s0 = growing.etat0[1]                                           # susceptible host plant density
    
    # solve the ODE problem for a first growing season
    problem = ODEProblem(growing.model, growing.etat0, tspan, growing.params, saveat=growing.pas)
    solution = solve(problem)
    # collect the results
    res.all_S = push!(res.all_S, solution[1, :])
    res.all_I = push!(res.all_I, solution[2, :])
    res.all_t = push!(res.all_t, solution.t)


    winter_length = growing.year - growing.tspan[2]

    # simulation for the rest of the time
    for _ in 1:years-1

        # update tspan of growing season 
        tspan = tspan .+ 365
        # collect the last values to get new initial conditions
        s_fin_g, i_fin_g = last(solution)
        # new initial conditions
        s0g = s0 * exp(-θ * π * exp(-μ * ((winter_length) * 1.0)) * i_fin_g / λ)
        i0g = s0 * (1 - exp(-θ * π * exp(-μ * ((winter_length) * 1.0)) * i_fin_g / λ))
        # encapsulation
        etat0 = @SVector [s0g, i0g]
        # solve the ODE problem for growing season
        problem = ODEProblem(growing.model, etat0, tspan, growing.params, saveat=growing.pas)
        solution = solve(problem)
        # collect the results
        res.all_S = push!(res.all_S, solution[1, :])
        res.all_I = push!(res.all_I, solution[2, :])
        res.all_t = push!(res.all_t, solution.t)

    end
    

    # convert days into years
    t = res.all_t ./ growing.year

    # plot S
    p1 = plot(t, res.all_S,
        label=false,
        xlims=[0, years],
        ylims=[0, s0],
        xlabel="Year",
        ylabel="\$S\$",
        c=:black)

    # plot I
    p2 = plot(t, res.all_I,
        label=false,
        xlims=[0, years],
        ylims=[0, s0 / 3],
        xlabel="Year",
        ylabel="\$I\$",
        c=:black)
    # plot S et I dans une même fenêtre
    plot(p1, p2,
        layout=(2, 1))
    title!("Simulation du modèle airborne compacte", subplot=1)

end



##############################################    TEST    ################################################################

# time
t_0 = 0
τ = 184
Τ = 365
t_transi = τ
t_fin = Τ
temps_simule = 5

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
paramsg = [α, β]
# others parameters
θ = 0.1
π = 1.0
μ = 0.0072
λ = 0.2938
others_params = [θ, π, μ, λ]

growing = Growing(params=paramsg, others_params=others_params, tspan=(t_0, t_transi))

simule(temps_simule, growing)