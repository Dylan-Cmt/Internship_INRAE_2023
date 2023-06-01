########################################################CO-oxydation (codim 2)###############################################################################
###################https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialCO/#CO-oxydation-(codim-2)###############################
#############################################################################################################################################################


using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# define the sup norm
norminf(x) = norm(x, Inf)

# vector field of the problem
function COm(z, p)
    @unpack r, K, c, h, b, m = p
    u1, u2 = z
    out = similar(z)
    out[1] = r * u1 * (1 - (u1 / K)) - (c * u1 * u2) / (h + u1)
    out[2] = b * u1 / (h + u1) * u2 - m * u2
    out
end

# parameters used in the model
par_com = (r=1.0, K=1.0, c=1.0, h=2.0, b=2.0, m=1.0)

recordCO(x, p) = (u2 = x[2])

# initial condition
z0 = [1.0, 2.5]

# Bifurcation Problem
prob = BifurcationProblem(COm, z0, par_com, (@lens _.K); recordFromSolution=recordCO)

# continuation parameters
opts_br = ContinuationPar(pMin=0.1, pMax=8.0, ds=0.01, dsmax=0.1)

# compute the branch of solutions
br = continuation(prob, PALC(), opts_br;
    plot=true, verbosity=2, normC=norminf)

# plot the branch
scene = plot(br, ylims=[-0.1, 5])