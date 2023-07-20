using Plots, DifferentialEquations, StaticArrays, Parameters, Test

# Setting up the problem

## TimeParam

"""
`TimeParam` stock toutes les variables relatives au temps
"""
@with_kw struct TimeParam
	T::Float64  = 365 ; @assert T > 0; @assert T <= 365 # year length
	τ::Float64  = 184 ; @assert τ <= T 					# crop season length
	Δt::Float64 = 0.1 ; @assert Δt > 0 					# step
	tspang::Tuple{Float64,Float64} = (0, τ) 			# tspan for growing
	tspanw::Tuple{Float64,Float64} = (τ, T) 			# tspan for winter	
end


## Param pour chaque modèle

"""
`Param` stock les paramètres biologiques, le type de modèle (compact/élaboré) et le nombre d'états du modèle.

Les différentes struct ont été trié à l'aide d'abstract types, de sorte à pouvoir utiliser des types plus généraux dans nos fonctions. 

Par exemple, `ParamAirborneElaborate1Strain` et `ParamSoilborneElaborate1Strain` ont été regroupé dans `Elaborate1Strain` car ils partagent la méthode `WinterSeason`.
"""
abstract type Param end


abstract type Elaborate1Strain <: Param end
abstract type Elaborate2Strains <: Param end
abstract type Compact1Strain <: Param end
abstract type Compact2Strains <: Param end


@with_kw struct ParamSoilborneCompact1Strain{T<:Float64} <: Compact1Strain
	# growing season parameters
	α::T  = 0.024   ; @assert α > 0 	# Infected host plants removal rate
	β::T  = 0.04875 ; @assert β > 0 	# Secondary infection rate

	# convertion into primary inoculum parameter
	θ::T  = 0.04875 ; @assert θ > 0 	# Primary infection rate
	Π::T  = 1.0 	; @assert Π > 0 	# Conversion rate from I to P (at the end of the season)
	μ::T  = 0.0072  ; @assert μ > 0
	λ::T  = 0.0072  ; @assert λ > 0
	n::T = 1.0 	  	; @assert n >= 0 	# Initial plant density

	# type of model and state number 
	isElaborate = false
	statelength = 2
end

@with_kw struct ParamAirborneElaborate1Strain{T<:Float64} <: Elaborate1Strain
	# growing season parameters
	Λ::T  = 0.052    ; @assert Λ > 0 	# Primary inoculum density independent depletion rate
	Θ::T  =0.04875  ; @assert Θ > 0 	# Primary infection rate
	α::T  = 0.024    ; @assert α > 0 	# Infected host plants removal rate
	β::T  = 0.04875 ; @assert β > 0 	# Secondary infection rate

	# convertion into primary inoculum parameter
	Π::T  = 1.0 	 ; @assert Π > 0 	# Conversion rate from I to P (at the end of the season)

	# winter-specific mortality parameter
	μ::T  = 0.0072  ; @assert μ > 0 	# Winter season mortality rate of primary inoculum

	# new susceptible host plant density parameter 
	n::T = 1.0 	  	 ; @assert n >= 0 	# Initial plant density

	# type of model and state number 
	isElaborate = true
	statelength = 3
end


## StateParam0

"""
`StateParam0` stock les états du modèle au début d'une saison.
"""
abstract type StateParam0 end

@with_kw struct StateCompact <: StateParam0
	S0::Float64 = 0.99 ; @assert S0 >= 0 
	I0::Float64 = 0.01 ; @assert I0 >= 0
	@assert S0+I0 <= 1
	State0 = @SVector [S0, I0] 					
end

@with_kw struct StateElaborate <: StateParam0
	P0::Float64 = 0.01 ; @assert P0 >= 0 
	S0::Float64 = 0.99 ; @assert S0 >= 0 
	I0::Float64 = 0.00 ; @assert I0 >= 0
	@assert S0+I0 <= 1
	State0 = @SVector [P0, S0, I0] 					
end


# Equations

"""
Chaque équation se différencie par le `Param` qu'il va prendre en argument.

Certaines équations peuvent prendre en argument un type plus large de `Param` car elles ont des paramètres en commun, comme par exemple le type `Elaborate1Strain` dans `WinterSeason`.
"""

## GrowingSeason
function GrowingSeason(State0::SVector,
						param::Compact1Strain,
						t::Real)
	S, I =  State0
	@unpack α, β = param
	
	dS = - β * S * I
	dI = + β * S * I - α * I
	
	@SVector [dS, dI]
end

function GrowingSeason(State0::SVector,
						param::ParamAirborneElaborate1Strain,
						t::Real)

	P, S, I = State0
	@unpack α, β, Λ, Θ = param
	
	dP = - Λ * P
	dS = - Θ * P * S - β * S * I
	dI = + Θ * P * S + β * S * I - α * I

	@SVector [dP, dS, dI]
end


## WinterSeason
function WinterSeason(State0::SVector,
					  param::Elaborate1Strain,
					  t::Real)
	P, S, I =  State0
	@unpack μ = param
	dP = −μ * P
	dS = 0
	dI = 0
	@SVector [dP, dS, dI]
end


# Seasons simulations

## Growing and Winter

"""
`growing` va prendre en argument un `StateParam0`, un `Param` et un `TimeParam`. Il va ensuite simuler une saison pour n'importe quel modèle et retourner une matrice contenant la simulation ainsi que les dernières valeurs de celle-ci.

`winter` va prendre en argument les dernières valeurs de la simulation de `growing`, un `Param` et un `TimeParam`. Il va ensuite simmuler l'hiver et retourner les mêmes objets que `growing`.

"""
function growing(sp::StateParam0,
				param::Param;
				tp::TimeParam=TimeParam())

	# simulation
	@unpack tspang, Δt = tp	
	prob = ODEProblem(GrowingSeason, sp.State0, tspang, param, saveat = Δt)
	sol  = solve(prob)

	# collect of last values
	res_end = last(sol)

	# build of results matrix
	res = []
	push!(res, sol.t)
	for i in 1:param.statelength
		# sol[i,:] = [...]
		push!(res, sol[i,:])
	end
	# res = [ [...], [...], ...]
	return res, res_end
end

function winter(res_end,
				param::Elaborate1Strain;
				tp::TimeParam=TimeParam())

	# compute new CI
	Pend, Send, Iend = res_end
	@unpack Π = param
	sp = StateElaborate(P0=Pend + Π*Iend, S0=0, I0=0)

	# simulation
	@unpack tspanw, Δt = tp
	prob = ODEProblem(WinterSeason, sp.State0, tspanw, param, saveat = Δt)
	sol  = solve(prob)

	# collect of last values
	res_end = last(sol)

	# build of results matrix
	res = []
	push!(res, sol.t)
	for i in 1:param.statelength
		# sol[i,:] = [...]
		push!(res, sol[i,:])
	end
	# res = [ [...], [...], ...]
	return res, res_end
end


## Years transitions

"""
D'une saison à l'autre, il est important de recalculer les nouvelles conditions initiales.

`yeartransition` va prendre en argument `res_end` (la fin d'une simulation), un `Param`, et éventuellement un `TimeParam`. Il retourne un nouvel état.
"""

# for a compact model
function yeartransition(res_end,
						param::ParamSoilborneCompact1Strain;
						tp::TimeParam=TimeParam())
	Send, Iend = res_end

	@unpack θ, Π, μ, λ, n = param
	@unpack T, τ = tp

	Snew = n * exp(-θ*Π*exp(-μ*(T-τ))/λ * Iend)
	Inew = n - Snew
	return StateCompact(S0=Snew, I0=Inew)
end

# for an elaborate model
function yeartransition(res_end,
						param::ParamAirborneElaborate1Strain;
						tp::TimeParam=TimeParam())
	Pend, Send, Iend = res_end
	@unpack n = param
	return StateElaborate(P0=Pend, S0=n ,I0=0)
end


# Problem solving during a year

"""
`simule` simule pendant 1 an le problème.
"""
function simule(sp::StateParam0,
				param::Param;
				tp::TimeParam=TimeParam())

	# simule growing and collect data as a vector of vectors
	res, res_end = growing(sp, param, tp=tp)	
	
	# if elaborate model: compute new CI and simule winter
	if param.isElaborate
		resw, res_end = winter(res_end, param, tp=tp)
		# add result to the growing simulation
		for i in eachindex(res)
			res[i] = vcat(res[i], resw[i])
		end
	end

	# compute new CI for growing season
	CI = yeartransition(res_end, param, tp=tp)
	
	return res, CI
end


# Problem solving during n years

"""
> `fill_mat` crée un matrice vide de `SVectors` avec des dimensions en adéquation avec la modélisation.
"""
function fill_mat(nyears::Int64,
					sp::StateParam0,
					param::Param;
					tp::TimeParam=TimeParam())
	
	@unpack T, τ, Δt = tp
	
	return Matrix{Vector{Float64}}(undef, nyears, length(sp.State0)+1)
end

"""
> `simule(nyears)` simule pendant `nyears` en appelant en boucle la fonction `simule`.
"""
function simule(nyears::Int64,
				sp::StateParam0,
				param::Param;
				tp::TimeParam=TimeParam())
	
	@test param.statelength==length(sp.State0)
	
	@unpack T, Δt = tp
	mat_res =  fill_mat(nyears, sp, param, tp=tp)

	CI = sp
	for i in 1:nyears
		res, CI = simule(CI, param, tp=tp)
		
		mat_res[i,:] = res
		mat_res[i,1] = mat_res[i,1] .+ (i-1)*T
	end
	
	return mat_res
end


# Plot


# parcourt un SVector et applique la transformation terme à terme
isWinter_vect(t,tp) =[mod(x, 1) < tp.τ/tp.T ? 0 : 1 for x in t]

# parcourt une ligne (vecteur de vecteurs) de la matrice
isWinter(t,tp) = [isWinter_vect(x,tp) for x in t[:,1]]


function Plots.plot(nyears::Int64,
					sp::StateCompact,
					param::Compact1Strain;
					tp::TimeParam=TimeParam())
	mat_res = simule(nyears, sp, param, tp=tp)

	# convert days into years
	t = mat_res[:,1] ./365 

	
    # plot S
    p1 = plot(t, mat_res[:,2],
        label=false,
        xlims=[0, nyears],
        ylims=[0, param.n],
        ylabel="\$S\$",
        c=:black)
	#p1 = plot!(t, isWinter(t,tp), fillrange = 0, fillcolor = :lightgray, fillalpha = 0.65, lw = 0, label=false, legend=:topright)

    # plot I
    p2 = plot(t, mat_res[:,3],
        label=false,
        xlims=[0, nyears],
        ylims=[0, param.n / 3],
        xlabel="Years",
        ylabel="\$I\$",
        c=:black)
	#p2 = plot!(t, isWinter(t,tp), fillrange = 0, fillcolor = :lightgray, fillalpha = 0.65, lw = 0, label=false, legend=:topright)
	
    # plot S et I dans une même fenêtre
    plot(p1, p2,
        layout=(2, 1))
    title!("Simulation du modèle airborne compacte", subplot=1)

end


function Plots.plot(nyears::Int64,
					sp::StateElaborate,
					param::Elaborate1Strain;
					tp::TimeParam=TimeParam())
	mat_res = simule(nyears, sp, param, tp=tp)

	# convert days into years
	t = mat_res[:,1] ./365 

	
    # plot S
    p1 = plot(t, mat_res[:,3],
        label=false, ylabel="\$S\$",
        xlims=[0, nyears], ylims=[0, param.n],
        c=:black, linestyle=:solid)
	#p1 = plot!(t, isWinter(t,tp), fillrange = 0, fillcolor = :lightgray, fillalpha = 0.65, lw = 0, label=false, legend=:topright)

    # plot I
    p2 = plot(t, mat_res[:,4],
        label=false, xlabel="Years", ylabel="\$I\$",
        xlims=[0, nyears], ylims=[0, param.n / 3],
        c=:black, linestyle=:solid)
	# plot P
    p2 = plot!(twinx(),t, mat_res[:,2],
        label=false, ylabel="\$P\$",
        xlims=[0, nyears], ylims=[0, param.n / 3],
		c=:black, linestyle=:dashdotdot)
	#p2 = plot!(t, isWinter(t,tp), fillrange = 0, fillcolor = :lightgray, fillalpha = 0.65, lw = 0, label=false, legend=:topright)
	
    # plot S et I dans une même fenêtre
    plot(p1, p2,
        layout=(2, 1))
    title!("Simulation du modèle soilborne élaboré", subplot=1)

end


# Test

# time parameter
tp = TimeParam()
# elaborate 1 strain model
paramE = ParamAirborneElaborate1Strain()
spE = StateElaborate()
# compact 1 strain model
paramC = ParamSoilborneCompact1Strain()
spC = StateCompact()


plot(4, spC, paramC)
plot(4, spE, paramE)
@time plot(100, spE, paramE);
