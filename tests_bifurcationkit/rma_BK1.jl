#####################################################pp2 example from AUTO07p (aBD + Hopf aBS)###############################################################
#############https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialPP2/#pp2-example-from-AUTO07p-(aBD-Hopf-aBS) ##################
#############################################################################################################################################################

using Revise, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)

# function to record information from a solution
recordFromSolution(x, p) = (u2 = x[2])

	
function pp2!(dz, z, p, t)
	@unpack r, K, c, h, b, m = p
	u1, u2 = z
	dz[1] = r*u1*(1- (u1/K)) - (c*u1*u2)/(h+u1) 
	dz[2] =	b*u1/(h+u1)*u2 - m*u2 
	dz
end
	
pp2(z, p) = pp2!(similar(z), z, p, 0)
	
# parameters of the model
par_pp2 = (r = 1.0, K = 1.0, c = 1.0, h = 2.0, b = 2.0, m = 1.0)

# initial condition
z0 = [1.0 , 2.5]
	
# bifurcation problem
prob = BifurcationProblem(pp2, z0, par_pp2,
	# specify the continuation parameter
	(@lens _.K), recordFromSolution = recordFromSolution)

# continuation options
opts_br = ContinuationPar(pMin = 0.1, pMax = 8.0, dsmax = 0.01,
	# options to detect bifurcations
	detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25,
	# number of eigenvalues
	nev = 2,
	# maximum number of continuation steps
	maxSteps = 1000,)

# parameter theta, see `? continuation`. Setting this to a non
# default value helps passing the transcritical bifurcation
diagram = bifurcationdiagram(prob, PALC(θ = 0.3),
	# very important parameter. It specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches of branches
	# at most in the present case.
	3,
	(args...) -> setproperties(opts_br; ds = 0.001, dsmax = 0.01, nInversion = 8, detectBifurcation = 3);
	# δp = -0.01,
	verbosity = 1, plot = true)
	

scene = plot(diagram; code = (), title="$(size(diagram)) branches", legend = true,
		ylims=[-0.1, 5])
