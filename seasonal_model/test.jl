# imports 
using DifferentialEquations                                         # for ODEProblem and solve
using Plots                                                         # for plot
using StaticArrays                                                  # for @SVector 
using Parameters                                                    # for @with_kw


##############################################    INITIALISATION    ################################################################

# time
t_0      = 0
τ        = 184                                                      # winter season duration
Τ        = 365                                                      # year duration
t_transi = Τ - τ                                                    # end of growing season / beginning of winter season
t_fin    = Τ                                                        # end of winter season / beginning of growing season 

#tspan
pas_t  = 1
tspang = (t_0, t_transi)
tspanw = (t_transi, t_fin)

# initial conditions
p0_g = 0.01                                                           # primary inoculum density
s0_g = 1.0                                                            # susceptible host plant density
i0_g = 0.0                                                            # infected host plant density
# encapsulation 
etat0 = @SVector [p0, s0, i0]

# accumulation of results
all_P, all_S, all_I = [], [], []


##############################################    SET UP FOR GROWING SEASON    ########################################################

function refreshInitialCond(mod::Growing)
    # collect the last data of winter season
    p_fin_w = last(all_P)                                          
    s_fin_w = last(all_S)
    i_fin_w = last(all_I)
    # initial conditions for growing season
    p_g = p_fin_w + π * i_fin_w                                    
    s_g = 0.0
    i_g = 0.0
    etat_g = @SVector [p_g, s_g, i_g]
    # refresh etat0 after the winter season
    mod.etat0 = etat_g
end

# parameters
α = 0.024                                                           # infected host plants removal rate per day
β = 0.04875                                                         # secondary infection rate per day per host plant unit
Λ = 0.052                                                           # within-season primary inoculum loss rate per day
Θ = 0.04875                                                         # primary infection rate per primary inoculum unit per day
paramsg = [α, β, Λ, Θ]

@with_kw struct Growing
    etat0::SVector{3,Float64} = @SVector [0.01, 1.0, 0.0]           # default value = initial value
    params::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    pas::Float64
    model::Function                                                 # ODE
end



# model for the growing season
function modelg(u::SVector{3,Float64}, params, t)
    α, β, Λ, Θ = params                                             # unpack the vectors into scalar
    p, s, i = u
    dp = -Λ * p                                                     # dot p
    ds = -Θ * p * s - β * s * i                                     # dot s
    di = Θ * p * s + β * s * i - α * i                              # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end


##############################################    SET UP FOR GROWING SEASON    ########################################################

π = 1                                                               # arbitrary primary inoculum unit per host plant unit
μ = 0.0072                                                          # per day

function refreshInitialCond(mod::Winter)
    # collect the last data of winter season
    p_fin_g = last(all_P)                                          
    s_fin_g = last(all_S)
    i_fin_g = last(all_I)
    # initial conditions for growing season
    p_w = p_fin_g                               
    s_w = s0_g
    i_w = 0.0
    etat_w = @SVector [p_w, s_w, i_w]
    # refresh etat0 after the winter season
    mod.etat0 = etat_w
end


@with_kw struct Winter
    etat0::SVector{3,Float64}
    params::Float64
    tspan::Tuple{Float64,Float64}
    pas::Float64
    model::Function                                                 # ODE
end

# model for the winter season
function modelw(u::SVector{3,Float64}, param, t)
    μ = param                                                      # unpack the vectors into scalar
    p, s, i = u
    dp = -μ * p                                                     # dot p
    ds = 0                                                          # dot s
    di = 0                                                          # dot i
    @SVector [dp, ds, di]                                           # return a new vector
end


##############################################    SOLVE AND RECORD OF SOLUTION    ########################################################

function record(mod::Union{Growing,Winter})
    problem = ODEProblem(mod.model, mod.etat0, mod.tspan, mod.params, saveat=mod.pas)
    solution = solve(problem)
    all_P = push!(all_P, solution[1, :])
    all_S = push!(all_S, solution[2, :])
    all_I = push!(all_I, solution[3, :])
end

