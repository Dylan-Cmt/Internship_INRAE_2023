#####################################################pp2 example from AUTO07p (aBD + Hopf aBS)###############################################################
#############https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialPP2/#pp2-example-from-AUTO07p-(aBD-Hopf-aBS) ##################
#############################################################################################################################################################

using Revise, Parameters, Setfield, Plots, LinearAlgebra# Setfield is for @lens
using BifurcationKit
const BK = BifurcationKit

# define the sup norm
norminf(x) = norm(x, Inf)

# function to record information from a solution
recordFromSolution(x, p) = (u2 = x[2])


function pp2!(dz, z, p, t)
    @unpack r, K, c, h, b, m = p
    u1, u2 = z
    dz[1] = r * u1 * (1 - (u1 / K)) - (c * u1 * u2) / (h + u1)
    dz[2] = b * u1 / (h + u1) * u2 - m * u2
    dz
end

pp2(z, p) = pp2!(similar(z), z, p, 0)

# parameters of the model
par_pp2 = (r=1.0, K=1.0, c=1.0, h=2.0, b=2.0, m=1.0)

# initial condition
z0 = [1.0, 2.5]

# bifurcation problem
prob = BifurcationProblem(pp2, z0, par_pp2,
    # specify the continuation parameter
    (@lens _.K), recordFromSolution=recordFromSolution)

# continuation options
opts_br = ContinuationPar(pMin=0.1, pMax=8.0, dsmax=0.01,
    # options to detect bifurcations
    detectBifurcation=3, nInversion=8, maxBisectionSteps=25,
    # number of eigenvalues
    nev=2,
    # maximum number of continuation steps
    maxSteps=1000,)

# parameter theta, see `? continuation`. Setting this to a non
# default value helps passing the transcritical bifurcation
diagram = bifurcationdiagram(prob, PALC(θ=0.3),
    # very important parameter. It specifies the maximum amount of recursion
    # when computing the bifurcation diagram. It means we allow computing branches of branches of branches
    # at most in the present case.
    3,
    (args...) -> setproperties(opts_br; ds=0.001, dsmax=0.01, nInversion=8, detectBifurcation=3);
    # δp = -0.01,
    verbosity=1, plot=true)


scene = plot(diagram; code=(), title="$(size(diagram)) branches", legend=true,
    ylims=[-0.1, 5])


# branch of the diagram with Hopf point
brH = BK.getBranch(diagram, (2, 1)).γ	# this line contains a BoundsError: attempt to access 0-element Vector{BifDiagNode} at index [1]

# newton parameters (with Default values)
optn_po = NewtonPar(tol=1e-8, maxIter=25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax=0.1, ds=0.01, dsmin=1e-4, newtonOptions=(@set optn_po.tol = 1e-8), tolStability=1e-2, detectBifurcation=1)

Mt = 101 # number of time sections
br_po = continuation(
    brH, 1, opts_po_cont,
    PeriodicOrbitTrapProblem(M=Mt;
        # specific linear solver for ODEs
        jacobian=:Dense);
    recordFromSolution=(x, p) -> (xtt = reshape(x[1:end-1], 2, Mt);
    return (max=maximum(xtt[1, :]),
        min=minimum(xtt[1, :]),
        period=x[end])),
    finaliseSolution=(z, tau, step, contResult; prob=nothing, kwargs...) -> begin
        # limit the period
        z.u[end] < 100
    end,
    normC=norminf)


plot(diagram, ylims=[-0.1, 5]);
plot!(br_po, label="Periodic orbits", legend=:bottomright);
