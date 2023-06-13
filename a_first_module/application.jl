include("DynamiquePopulation.jl")
# To load a module from a locally defined module, a dot needs to be added before the module name
using .DynamiquePopulation                      # import .DynamiquePopulation as DP OK mais pas DynamiquePopulation as DP
using StaticArrays

function main()

    # tspan
    t_0 = 0
    t_fin = 30.0
    pas = 0.01
    tspan = (t_0, t_fin)

    # parameters
    r = 1.0
    c = 1.0
    b = 1.0
    m = 1.0
    params_LV = [r, c, b, m]

    # initial density of population
    #etat01 = @SVector [0.1]
    etat02 = @SVector [1.0, 2.5]

    # modeling
    #=
    # Populations isolées
    mod_malthus = DynamiquePopulation.Modeling(etat01, [3, 2], (0.0, 20.0), 0.01, malthus)
    mod_logistic = DynamiquePopulation.Modeling(etat01, [1, 10], (0.0, 10.0), 0.01, logistic)
    mod_allee = DynamiquePopulation.Modeling(etat01, [1, 10, 2], (0.0, 3.0), 0.01, allee)
    
    # Populations en intéraction
    mod_rma = DynamiquePopulation.Modeling(etat02, [1.0, 10.0, 1.0, 2.0, 2.0, 1.0], (0.0, 80.0), 0.01, rma)
    =#
    mod_lv = DynamiquePopulation.Modeling(etat02, params_LV, tspan, pas, DynamiquePopulation.lv)
    
    # test
    DynamiquePopulation.simule(mod_lv,
        plot=true,
        xlabel="Temps \$t\$",
        ylabel="Densités de population",
        label=["\$x(t)\$" "\$y(t)\$"],
        title="Simulation du modèle de Lotka Volterra ")

    #DynamiquePopulation.simule(mod_lv)

end

main()