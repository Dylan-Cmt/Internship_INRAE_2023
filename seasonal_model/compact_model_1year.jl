# imports 
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 

# time
t_0      = 0
τ        = 184                                                             
Τ        = 365                                                            
t_transi = Τ - τ
t_fin    = Τ

#tspan
pas_t  = 1
tspang = (t_0, t_transi)
tspanw = (t_transi, t_fin)

##############################################    GROWING SEASON: year 1    ################################################################


# initial conditions
s0 = 1.0                                                    # susceptible host plant density
i0 = 0.0                                                    # infected host plant density
# encapsulation 
etat0 = @SVector [s0, i0]

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
paramsg = [α, β]


# model for the growing season
function modelg(u::SVector{2,Float64}, params, t)
    α, β = params                                                   # unpack the vectors into scalar
    s, i = u
    ds   = - β * s * i                                              # dot s
    di   = + β * s * i - α * i                                      # dot i
    @SVector [ds, di]                                               # return a new vector
end

problemg  = ODEProblem(modelg, etat0, tspang, paramsg, saveat=pas_t)
solutiong = solve(problemg)





##############################################    WINTER SEASON: year 1    ################################################################
 
# nothing is happening during winter ...

##############################################    GROWING SEASON: year 2    ################################################################


# winter season data recovery
s_fin_g, i_fin_g = last(solutiong)

# new initial conditions
s1 = s0 * exp(-θ * π * exp(-μ(Τ - τ)) * i_fin_g / λ)
i1 = s0*(1 - exp(-θ * π * exp(-μ(Τ - τ)) * i_fin_g / λ))
# encapsulation 
etat0w = @SVector [s1, i1]

problemg = ODEProblem(modelg, etat1, tspang .+ 360, paramsg, saveat=pas_t)
solutiong = solve(problemg)

##############################################    GROWING SEASON: year 2    ################################################################

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