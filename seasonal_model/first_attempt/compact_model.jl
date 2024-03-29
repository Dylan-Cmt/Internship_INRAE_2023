#################################################     IMPORTS    ###########################################################################
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 
using Parameters                                                    # for @with_kw
###############################################    STRUCTS AND FUNCTIONS    ################################################################

"""
    struct Growing

Contains the model informations for the growing season,
with few default values.
"""
@with_kw struct Growing
    etat0::SVector{2,Float64}
    params::Vector{Float64}
    tspan::Tuple{Int64,Int64}
    pas = 1
    model::Function = modelg
end

"""
    struct OtherParameters

Contains other parameters.
"""
mutable struct OtherParameters
    params::Union{Int64,Float64,Vector{Float64}}
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
function simule(years, growing::Growing, other::OtherParameters)

    # Creat a Result type to collect results
    res = Result()

    # collect tspan
    tspan = growing.tspan

    # collect others parameters
    θ, π, μ, λ = other.params

    # initial condition
    s0 = growing.etat0[1]                                           # susceptible host plant density
    
    # solve the ODE problem for a first growing season
    problem = ODEProblem(growing.model, growing.etat0, tspan, growing.params, saveat=growing.pas)
    solution = solve(problem)
    # collect the results
    res.all_S = push!(res.all_S, solution[1, :])
    res.all_I = push!(res.all_I, solution[2, :])
    res.all_t = push!(res.all_t, solution.t)

    Τ = 365
    τ = growing.tspan[2]

    # make a first strip
    v1, v2 = [[0, τ]], [[0, τ]]
   
    # simulation for the rest of the time
    for i in 1:1:years-1

        # complete the vector to make the others stips
        u = [i * Τ, i * Τ + τ]
        push!(v1, u)
        push!(v2, u)

        # update tspan of growing season 
        tspan = tspan .+ 365
        # collect the last values to get new initial conditions
        s_fin_g, i_fin_g = last(solution)
        # new initial conditions
        s0g = s0 * exp(-θ * π * exp(-μ * ((Τ - τ) * 1.0)) * i_fin_g / λ)
        i0g = s0 * (1 - exp(-θ * π * exp(-μ * ((Τ - τ) * 1.0)) * i_fin_g / λ))
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
    t = res.all_t ./ Τ
    v1 = v1 ./ Τ
    v2 = v2 ./ Τ

    # add stips to plot p1
    p1 = vspan(v1[1],
        color=:lightgray,
        label="growing season")
    p1 = vspan!(v1[2:end],
        color=:lightgray,
        label=false)
    # plot S
    p1 = plot!(t, res.all_S,
        label=false,
        xlims=[0, years],
        ylims=[0, s0],
        ylabel="\$S\$",
        c=:black)

    # add stips to plot p2
    p2 = vspan(v2[1],
        color=:lightgray,
        label="growing season")
    p2 = vspan!(v2[2:end],
        color=:lightgray,
        label=false)
    # plot I
    p2 = plot!(t, res.all_I,
        label=false,
        xlims=[0, years],
        ylims=[0, s0 / 3],
        xlabel="Years",
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
tspang = (t_0, t_transi)
temps_simule = 5

# initial conditions
etat0 = @SVector [1.0, 0.01]

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

growing = Growing(etat0=etat0, params=paramsg, tspan=tspang)
other = OtherParameters(others_params)

simule(temps_simule, growing, other)