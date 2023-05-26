#####################################################Neural mass equation (Hopf aBS)#########################################################################
################https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE/#Neural-mass-equation-(Hopf-aBS)#######################
#############################################################################################################################################################

using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

# vector field
function TMvf!(dz, z, p, t)
    @unpack r, K, c, h, b, m = p
    u1, u2 = z
    dz[1] = r * u1 * (1 - (u1 / K)) - (c * u1 * u2) / (h + u1)
    dz[2] = b * u1 / (h + u1) * u2 - m * u2
    dz
end

# out of place method
TMvf(z, p) = TMvf!(similar(z), z, p, 0)

# parameter values
par_tm = (r=1.0, K=1.0, c=1.0, h=2.0, b=2.0, m=1.0)

# initial condition
z0 = [1.0, 2.5]

# Bifurcation Problem
prob = BifurcationProblem(TMvf, z0, par_tm, (@lens _.K);
    recordFromSolution=(x, p) -> (u2 = x[2]))


# continuation options
opts_br = ContinuationPar(pMin=.1, pMax=8.0,
    # parameters to have a smooth result
    ds=0.009, dsmax=0.01,)

# continuation of equilibria
br = continuation(prob, PALC(θ=0.3), opts_br;
    plot=true, normC=norminf)

scene = plot(br, plotfold=false, markersize=3, legend=:topleft)