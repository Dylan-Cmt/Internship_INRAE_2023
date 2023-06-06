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


# accumulation of results (also type missing to show a discontinuity)
all_P::Vector{Union{Missing,Float64}} = []
all_S::Vector{Union{Missing,Float64}} = []
all_I::Vector{Union{Missing,Float64}} = []
all_t::Vector{Union{Missing,Float64}} = []


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

# put solution's values somewhere to plot later
all_P = vcat(all_P, solution[1, :])
all_S = vcat(all_S, solution[2, :])
all_I = vcat(all_I, solution[3, :])
all_t = vcat(all_t, solution.t)


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

# put solution's values somewhere to plot later
all_P = vcat(all_P, solutionw[1, :])
all_S = vcat(all_S, solutionw[2, :])
all_I = vcat(all_I, solutionw[3, :])
all_t = vcat(all_t, solutionw.t)
##############################################    GROWING SEASON: year 2    ################################################################

#=

p1 = plot(all_t, all_I,
    label=false,
    c=:red,
    xlabel="Days",
    ylabel="\$I(t)\$",
    linestyle=:solid,
    xlims=[0, 365],
    ylims=[0, 1 / 3])

p1 = plot!(twinx(), all_t, all_P,
    c=:black,
    label=false,
    ylabel="\$P(t)\$",
    size=(400, 300),
    linestyle=:dashdotdot,
    ylims=[0, 1 / 3])

p2 = plot(all_t, all_S, xlims=[0, 365], ylims=[0, s0], label=false, ylabel="\$S(t)\$")
plot(p2, p1, layout=(2,1))




# collect winter season data
p_fin_w, s_fin_w, i_fin_w = last(solutionw)

# new initial conditions
p0g = p_fin_w
s0g = s0                                                          
i0g = 0.0
# encapsulation 
etat0w = @SVector [p0g, s0g, i0g]
=#