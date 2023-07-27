using Parameters, StaticArrays, AxisArrays, Plots, DifferentialEquations, Test

"""
    GrowingSeason(State0::SVector, param::Compact1Strain, t::Real)

This is the function to enter in ODEProblem from DifferentialEquations.jl. 
	
For a compact model, it returns the ODE associated to the growing season.

See also [`WinterSeason`](@ref), [`ODEProblem`](@ref), [`growing`](@ref).
"""
function GrowingSeason(State0::SVector,
    param::Compact1Strain,
    t::Real)
    S, I = State0
    @unpack α, β = param

    dS = -β * S * I
    dI = +β * S * I - α * I

    @SVector [dS, dI]
end

"""
    GrowingSeason(State0::SVector, param::Compact1Strain, t::Real)

This is the function to enter in ODEProblem from DifferentialEquations.jl. 
	
For an elaborate model, it returns the ODE associated to the growing season.

See also [`WinterSeason`](@ref), [`ODEProblem`](@ref), [`growing`](@ref).
"""
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

"""
    WinterSeason(State0::SVector, param::Elaborate1Strain, t::Real)

This is the function to enter in ODEProblem from DifferentialEquations.jl. 
	
It only exists for elaborate models, and it returns the ODE associated to the winter season.

See also [`GrowingSeason`](@ref), [`ODEProblem`](@ref), [`winter`](@ref).
"""
function WinterSeason(State0::SVector,
					  param::Elaborate1Strain,
					  t::Real)
	P, S, I =  State0
	@unpack μ = param
	dP = -μ * P
	dS = 0
	dI = 0
	@SVector [dP, dS, dI]
end

"""
    growing(sp::StateParam0, param::Param; tp::TimeParam=TimeParam())

Simulates the growing season for any model, using ODEProblem from DifferentialEquations.jl.

It returns a vector of vectors that contains the simulation for a season, and also the last values of the simulation.

See also [`GrowingSeason`](@ref), [`ODEProblem`](@ref), [`winter`](@ref), [`simule`](@ref).
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

"""
    winter(res_end,	param::Elaborate1Strain; tp::TimeParam=TimeParam())

Compute new initial conditions from the last values of the previous growing season.
Then it simulates the winter season for elaborate 1 strain models, using ODEProblem from DifferentialEquations.jl.

It returns a vector of vectors that contains the simulation, and the last values of the simulation.

See also [`WinterSeason`](@ref), [`ODEProblem`](@ref), [`winter`](@ref), [`simule`](@ref).
"""
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

"""
    yeartransition(res_end,	param::ParamSoilborneCompact1Strain; tp::TimeParam=TimeParam())

Compute new initial conditions from the last values of growing season simulation.

It returns a StateCompact object.

See also [`simule`](@ref).
"""
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

"""
    yeartransition(res_end,	param::ParamAirborneElaborate1Strain; tp::TimeParam=TimeParam())

Compute new initial conditions from the last values of winter season simulation.

It returns a StateElaborate object.

See also [`simule`](@ref).
"""
function yeartransition(res_end,
						param::ParamAirborneElaborate1Strain;
						tp::TimeParam=TimeParam())
	Pend, Send, Iend = res_end
	@unpack n = param
	return StateElaborate(P0=Pend, S0=n ,I0=0)
end

"""
    simule(sp::StateParam0,	param::Param; tp::TimeParam=TimeParam())

Simulate a year for any model.

It returns a vector of vectors that contains one year of simulation, and also the last values of the simulation.

See also [`growing`](@ref), [`winter`](@ref), [`yeartransition`](@ref), [`fill_mat`](@ref), [`affiche`](@ref).
"""
function simule(sp::StateParam0,
				param::Param;
				tp::TimeParam=TimeParam())

	# simulate growing and collect data as a vector of vectors
	res, res_end = growing(sp, param, tp=tp)	
	
	# if elaborate model: compute new CI and simulate winter
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

"""
    fill_mat(nyears::Int64,	sp::StateParam0, param::Param; tp::TimeParam=TimeParam())

Construct an empty named matrix to stock nyears of simulation.

Labels are filled automatically.

See also [`simule`](@ref).
"""
function fill_mat(nyears::Int64,
					sp::StateParam0,
					param::Param;
					tp::TimeParam=TimeParam())
	
	@unpack T, τ, Δt = tp

	# autofill axis
	years = Symbol.(["annee$i" for i in 1:nyears])
	col = [:time]
	for i in 1:length(fieldnames(typeof(sp)))-1
		push!(col, fieldnames(typeof(sp))[i])
	end
	
	# creat undef matrix
	mat = Matrix{Vector{Float64}}(undef, nyears, length(sp.State0)+1)
	
	return AxisArray(mat, Axis{:lignes}(years), Axis{:colonnes}(col))
end

"""
    simule(nyears::Int64, sp::StateParam0, param::Param; tp::TimeParam=TimeParam())

Simulate n years for any model.

It returns a named matrix that contains n years of simulation.
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

"""
    isWinter(t,tp)

From a time vector, returns a new vector of 0 and 1 to identify growing and winter seasons.
"""
isWinter(t, tp) = [mod(x, 1) < tp.τ / tp.T ? 0 : 1 for x in t]

"""
    affiche(nyears::Int64, sp::StateParam0, param::Param; tp::TimeParam=TimeParam())

Make a simulation of n years for any model.
Plot the solutions of this simulation.

See also [`simule`](@ref).
"""
function affiche(nyears::Int64,
    sp::StateParam0,
    param::Param;
    tp::TimeParam=TimeParam())
    # simulate
    mat = simule(nyears, sp, param)

	t = mat[:, :time] ./ 365
    simulTime = 0:tp.Δt*0.001:nyears

    # plot S0
    p1 = plot(t, mat[:, :S0], label=false, ylabel="S", c=:black, linestyle=:solid)
    # add stripes
    p1 = plot!(simulTime, isWinter(simulTime, tp), fillrange=0, fillcolor=:lightgray, fillalpha=0.65, lw=0, label="winter")

    # plot everything else
    p2 = plot()
    for i in 2:size(mat)[2]
        if mat[:, i] != mat[:, :S0]
            p2 = plot!(t, mat[:, i]
                #, label = String(fieldnames(typeof(sp))[i-1]))
                , label=false, ylims=[0, param.n / 3], c=:black, linestyle=:solid)
        end
    end
    # add stripes
    p2 = plot!(simulTime, isWinter(simulTime, tp), fillrange=0, fillcolor=:lightgray, fillalpha=0.65, lw=0, label="winter", xlabel="Years")
	
    # plot S and everything else in two subplots
    plot(p1, p2,
        layout=(2, 1))
end