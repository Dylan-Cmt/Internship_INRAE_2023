################################################Extended Lorenz-84 model (codim 2 + BT/ZH aBS)###############################################################
############https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS)##################
#############################################################################################################################################################

using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

# vector field
function Lor(z, p)
    @unpack r, K, c, h, b, m = p
    u1, u2 = z
    [
        r*u1*(1- (u1/K)) - (c*u1*u2)/(h+u1),
        b * u1 / (h + u1) * u2 - m * u2
    ]
end

# parameter values
parlor = (r=1.0, K=1.0, c=1.0, h=2.0, b=2.0, m=1.0)

# initial condition
z0 = [1.0, 2.5]

# bifurcation problem
recordFromSolutionLor(x, p) = (u2 = x[2])

prob = BifurcationProblem(Lor, z0, setproperties(parlor; K=8.0), (@lens _.K);
    recordFromSolution=recordFromSolutionLor)

# continuation options
opts_br = ContinuationPar(pMin=0.1, pMax=8.0, ds=0.01, dsmax=0.1,
    # Optional: bisection options for locating bifurcations
    nInversion=8, maxBisectionSteps=25,
    # number of eigenvalues
    nev=2, maxSteps=1000)

# compute the branch of solutions
br = @time continuation(prob, PALC(), opts_br;
    normC=norminf,
    bothside=true)

scene = plot(br, plotfold=false, markersize=4, legend=:topleft, ylims=[-0.1, 5])