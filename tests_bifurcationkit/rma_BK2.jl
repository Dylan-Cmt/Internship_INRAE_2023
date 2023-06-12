#####################################################Neural mass equation (Hopf aBS)#########################################################################
################https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE/#Neural-mass-equation-(Hopf-aBS)#######################
#############################################################################################################################################################

using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

########################################################## encode the ODE ######################################################

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
par_tm = (r=1.0, K=8.0, c=1.0, h=2.0, b=2.0, m=1.0)

# initial condition
z0 = [1.0, 2.5]

# Bifurcation Problem
prob = BifurcationProblem(TMvf, z0, par_tm, (@lens _.K);
    recordFromSolution=(x, p) -> (u2 = x[2]))

########################################################### compute the branch of equilibria ################################

# continuation options
opts_br = ContinuationPar(pMin=0.1, pMax=8.0,
    # parameters to have a smooth result
    ds=-0.04, dsmax=0.1,)

# continuation of equilibria
br = continuation(prob, PALC(θ=0.3), opts_br;
    plot=true, normC=norminf)

scene = plot(br, plotfold=false, markersize=3, legend=:topleft, ylims=[-0.1, 5])

################################################## Branch of periodic orbits with Trapezoid method #########################

# newton parameters
optn_po = NewtonPar(verbose=true, tol=1e-10, maxIter=10)

# continuation parameters
opts_po_cont = ContinuationPar(opts_br, dsmax=0.1, ds=-0.0001, dsmin=1e-4,
    maxSteps=90, newtonOptions=(@set optn_po.tol = 1e-7), tolStability=1e-8)

# arguments for periodic orbits
args_po = (recordFromSolution=(x, p) -> begin
        xtt = BK.getPeriodicOrbit(p.prob, x, p.p)
        return (max=maximum(xtt[2, :]),
            min=minimum(xtt[2, :]),
            period=getPeriod(p.prob, x, p.p))
    end,
    plotSolution=(x, p; k...) -> begin
        xtt = BK.getPeriodicOrbit(p.prob, x, p.p)
        arg = (marker=:d, markersize=1)
        plot!(xtt.t, xtt[1, :]; label="x", arg..., k...)
        plot!(xtt.t, xtt[2, :]; label="y", arg..., k...)
        plot!(br; subplot=1, putspecialptlegend=false)
    end,
    # we use the supremum norm
    normC=norminf)

Mt = 200 # number of time sections
br_potrap = continuation(
    # we want to branch form the 1th bif. point
    br, 1, opts_po_cont,
    # we want to use the Trapeze method to locate PO
    PeriodicOrbitTrapProblem(M=Mt);
    # regular continuation options
    verbosity=2, plot=true,
    args_po...
)

scene = plot(br, br_potrap, markersize=3)
plot!(scene, br_potrap.param, br_potrap.min, label="")
ylims!(0, 8)