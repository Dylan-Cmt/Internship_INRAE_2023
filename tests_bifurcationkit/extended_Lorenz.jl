using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)


################################## Problem setting ####################################################

# vector field
function Lor(u, p)
    @unpack α, β, γ, δ, G, F, T = p
    X, Y, Z, U = u
    [
        -Y^2 - Z^2 - α * X + α * F - γ * U^2,
        X * Y - β * X * Z - Y + G,
        β * X * Y + X * Z - Z,
        -δ * U + γ * U * X + T
    ]
end

# parameter values
parlor = (α=1 // 4, β=1, G=0.25, δ=1.04, γ=0.987, F=1.7620532879639, T=0.0001265)

# initial condition
z0 = [2.9787004394953343, -0.03868302503393752, 0.058232737694740085, -0.02105288273117459]



################################## Continuation and codim 1 bifurcations ##############################

# bifurcation problem
recordFromSolutionLor(x, p) = (X=x[1], Y=x[2], Z=x[3], U=x[4])
prob = BifurcationProblem(Lor, z0, setproperties(parlor; T=0.04, F=3.0), (@lens _.F);
    recordFromSolution=recordFromSolutionLor)

# continuation options
opts_br = ContinuationPar(pMin=-1.5, pMax=3.0, ds=0.002, dsmax=0.15,
    # Optional: bisection options for locating bifurcations
    nInversion=6, maxBisectionSteps=25,
    # number of eigenvalues
    nev=4, maxSteps=200)

# compute the branch of solutions
br = @time continuation(prob, PALC(), opts_br;
    normC=norminf,
    bothside=true)

scene = plot(br, plotfold=false, markersize=4, legend=:topleft)

################################## Continuation of Fold points ########################################

# function to record the current state
sn_codim2 = continuation(br, 5, (@lens _.T), ContinuationPar(opts_br, pMax=3.2, pMin=-0.1, detectBifurcation=1, dsmin=1e-5, ds=-0.001, dsmax=0.005, nInversion=10, maxSteps=130, maxBisectionSteps=55); normC=norminf,
    # detection of codim 2 bifurcations with bisection
    detectCodim2Bifurcation=2,
    # we update the Fold problem at every continuation step
    updateMinAugEveryStep=1,
    startWithEigen=false,
    # we save the different components for plotting
    recordFromSolution=recordFromSolutionLor
)

scene = plot(sn_codim2, vars=(:X, :U), branchlabel="Folds", ylims=(-0.5, 0.5))

################################## Continuation of Hopf points ########################################

hp_codim2_1 = continuation((@set br.alg.tangent = Bordered()), 3, (@lens _.T), ContinuationPar(opts_br, ds=-0.001, dsmax=0.02, dsmin=1e-4, nInversion=6, detectBifurcation=1); normC=norminf,
    # detection of codim 2 bifurcations with bisection
    detectCodim2Bifurcation=2,
    # we update the Fold problem at every continuation step
    updateMinAugEveryStep=1,
    # we save the different components for plotting
    recordFromSolution=recordFromSolutionLor,
    # compute both sides of the initial condition
    bothside=true
)

plot(sn_codim2, vars=(:X, :U), branchlabel="Folds", xlims=[1, 1.3], ylims=[-0.7, 0.7])
plot!(hp_codim2_1, vars=(:X, :U), branchlabel="Hopfs")

########################## Continuation of Hopf points from the Bogdanov-Takens point ##################

hp_from_bt = continuation((@set sn_codim2.alg.tangent = Bordered()), 4, ContinuationPar(opts_br, ds=-0.001, dsmax=0.02, dsmin=1e-4,
        nInversion=6, detectBifurcation=1); normC=norminf,
    # detection of codim 2 bifurcations with bisection
    detectCodim2Bifurcation=2,
    # we update the Fold problem at every continuation step
    updateMinAugEveryStep=1,
    # we save the different components for plotting
    recordFromSolution=recordFromSolutionLor
)

plot(sn_codim2, vars=(:X, :U), branchlabel="SN", xlims=[0.95, 1.3], ylims=[-0.7, 0.7])
plot!(hp_codim2_1, vars=(:X, :U), branchlabel="Hopf1")
plot!(hp_from_bt, vars=(:X, :U), branchlabel="Hopf2")

########################## Continuation of Hopf points from the Zero-Hopf point ########################

hp_from_zh = continuation((@set sn_codim2.alg.tangent = Bordered()), 2, ContinuationPar(opts_br, ds=0.001, dsmax=0.02, dsmin=1e-4, nInversion=6, detectBifurcation=1, maxSteps=150);
    normC=norminf,
    detectCodim2Bifurcation=2,
    updateMinAugEveryStep=1,
    startWithEigen=true,
    recordFromSolution=recordFromSolutionLor,
    bothside=false,
    bdlinsolver=MatrixBLS()
)

plot(sn_codim2, vars=(:X, :U), xlims=[0.95, 1.3], ylims=[-0.7, 0.75])
plot!(hp_codim2_1, vars=(:X, :U), branchlabel="Hopf")
plot!(hp_from_bt, vars=(:X, :U), branchlabel="Hopf2")
plot!(hp_from_zh, vars=(:X, :U), branchlabel="Hopf", plotspecialpoints=false, legend=:topleft)