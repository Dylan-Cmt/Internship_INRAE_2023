# imports 
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 

# time
t_0      = 0
τ        = 184                                                      # days
Τ        = 365                                                      # days
t_transi = Τ - τ
t_fin    = Τ

#tspan
pas_t  = 1
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
params = [α, β, Λ, Θ]


# model for the growing season
function modelg(u::SVector{3,Float64}, params, t)
    α, β, Λ, Θ = params                                             # unpack the vectors into scalar
    p, s, i    = u
    dp         = -Λ * p                                             # dot p
    ds         = -Θ * p * s - β * s * i                             # dot s
    di         =  Θ * p * s + β * s * i - α * i                     # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

problem  = ODEProblem(modelg, etat0, tspang, params, saveat=pas_t)
solution = solve(problem)
# plot S
p1 = plot(solution.t, [v[2] for v in solution.u],label=false,
    xlims  = [0, Τ],
    ylims  = [0, s0+0.2],
    xlabel = "Year",
    ylabel = "\$S\$")
# plot I
p2 = plot(solution.t, [v[3] for v in solution.u],label=false,
    xlims  = [0, Τ],
    ylims  = [0, s0/3],
    xlabel = "Year",
    ylabel = "\$I\$")
# plot S et I dans une même fenêtre
plot(p1, p2,
    layout = (2, 1))
title!("Simulation du modèle airborne élaboré",subplot=1)


##############################################    WINTER SEASON: year 1    ################################################################



# collect growing season data
p_fin_g, s_fin_g, i_fin_g = last(solution)

# additional parameter
π = 1                                                               # arbitrary primary inoculum unit per host plant unit
μ = 0.0072                                                          # per day


# new initial conditions
p0w = p_fin_g + π * i_fin_g
s0w = 0.0                                                           # arbitrary host plant unit
i0w = 0.0
# encapsulation 
etat0w = @SVector [p0w, s0w, i0w]

# model for the winter season
function modelw(u::SVector{3,Float64}, params, t)
    μ       = params                                                # unpack the vectors into scalar
    p, s, i = u
    dp      = -μ * p                                                # dot p
    ds      = 0                                                     # dot s
    di      = 0                                                     # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

problemw  = ODEProblem(modelw, etat0w, tspanw, μ, saveat=pas_t)
solutionw = solve(problemw)


# plot S (the only one which is non-zero)
plot!(solutionw,
    xlims  = [0, Τ],
    ylims  = [0, s0 + 0.2],
    xlabel = "Year",
    ylabel = "\$S\$")
title!("Simulation du modèle airborne élaboré", subplot=1)
# plot I
p2 = plot(solutionw.t, zeros(length(solutionw.t)) .+ 1, label=false,
    xlims  = [0, Τ],
    ylims  = [0, s0/3],
    xlabel = "Year",
    ylabel = "\$I\$",
    c=:red)
# plot S et I dans une même fenêtre
plot(p1, p2,
    layout = (2, 1))
title!("Simulation du modèle airborne élaboré",subplot=1)

##############################################    GROWING SEASON: year 2    ################################################################


# collect winter season data
p_fin_w, s_fin_w, i_fin_w = last(solutionw)

# new initial conditions
p0g = p_fin_w
s0g = s0                                                          
i0g = 0.0
# encapsulation 
etat0w = @SVector [p0g, s0g, i0g]



"""
@with_kw struct Mod
    etat0::Union{SVector{1,Float64},SVector{2,Float64}}
    p::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas::Float64
    f::Function                                                     # modelg or modelw
    season::Bool                                                    # true or false for growing or winter respectively
end

simule(years)                                                       # simulates (and plot?) a succession of growing and winter season
    loop
        # growing
        setInitialCond(m::Mod)                                      # collects p,s and i at the end of the season and returns the initial conditions encapsulated
        modelg(m::Mod)                                              # set the problem
        solve(m::Mod)                                               # uses ODEProblem then solve, put plots in variables

        #winter
        Mod.season = false                                          # changes the season
        setInitialCond(m::Mod)
        modelg(m::Mod)
        solve(m::Mod)
    
    display all plots 
"""