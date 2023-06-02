# imports 
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 

# time
t_0 = 0
τ   = 184                                                           # days
Τ = 365                                                             # days
t_transi = Τ - τ
t_fin = Τ

#tspan
pas_t = 1
tspang = (t_0, t_transi)
tspanw = (t_transi, t_fin)

##############################################    GROWING SEASON: year 1    ################################################################


# initial conditions
p0g = 0.01                                                          # primary inoculum density
s0g = 1.0                                                           # susceptible host plant density
i0g = 0.0                                                           # infected host plant density
# encapsulation 
etat0g = @SVector [p0g, s0g, i0g]

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
Λ = 0.052                                                           # within-season primary inoculum loss rate per day
Θ = 0.04875                                                         # primary infection rate per primary inoculum unit per day
paramsg = [α, β, Λ, Θ]


# model for the growing season
function modelg(u::SVector{3,Float64}, params, t)
    α, β, Λ, Θ = params                                             # unpack the vectors into scalar
    p, s, i    = u
    dp         = -Λ * p                                             # dot p
    ds         = -Θ * p * s - β * s * i                             # dot s
    di         =  Θ * p * s + β * s * i - α * i                     # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

problemg  = ODEProblem(modelg, etat0g, tspang, paramsg, saveat=pas_t)
solutiong = solve(problemg)
# plot S
p1 = plot(solutiong.t, [v[2] for v in solutiong.u],label=false,
    xlims=[0, Τ],
    ylims=[0, s0g+0.2],
    xlabel="Year",
    ylabel="\$S\$")
# plot I
p2 = plot(solutiong.t, [v[3] for v in solutiong.u],label=false,
    xlims=[0, Τ],
    ylims=[0, s0g/3],
    xlabel="Year",
    ylabel="\$I\$")
# plot S et I dans une même fenêtre
plot(p1, p2,
    layout=(2, 1))
title!("Simulation du modèle airborne élaboré",subplot=1)



##############################################    WINTER SEASON: year 1    ################################################################



# growing season data recovery
p_fin_g, s_fin_g, i_fin_g = last(solutiong)

# additional parameter
π = 1                                                               # arbitrary primary inoculum unit per host plant unit
μ = 0.0072                                                          # per day


# new initial conditions
p0w = p_fin_g + π * i_fin_g
s0w = 0.0                                                           # arbitrary host plant unit
i0w = 0.0
# encapsulation 
etat0w = @SVector [p0w, s0w, i0w]

# model for the growing season
function modelw(u::SVector{3,Float64}, params, t)
    μ       = params                                                # unpack the vectors into scalar
    p, s, i = u
    dp      = -μ * p                                                # dot p
    ds      = 0                                                     # dot s
    di      = 0                                                     # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end

problemw = ODEProblem(modelw, etat0w, tspanw, μ, saveat=pas_t)
solutionw = solve(problemw)


# plot S (the only one which is non-zero)
plot(solutionw,
    xlims=[0, Τ],
    ylims=[0, s0g + 0.2],
    xlabel="Year",
    ylabel="\$S\$")
title!("Simulation du modèle airborne élaboré", subplot=1)

##############################################    GROWING SEASON: year 2    ################################################################


# winter season data recovery
p_fin_w, s_fin_w, i_fin_w = last(solutionw)

# new initial conditions
p1g = p_fin_w + π * i_fin_w
s1g = 0.0                                                          
i1g = 0.0
# encapsulation 
etat0w = @SVector [p1g, s1g, i1g]