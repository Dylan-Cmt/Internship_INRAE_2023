# To load a module from a locally defined module, a dot needs to be added before the module name
using .DynamiquePopulation                  # import .DynamiquePopulation as DP OK mais pas DynamiquePopulation as DP
using .Modeles                                 # toutes les fonctions de EDOs ont été exportées
using StaticArrays

function main()

    # initial density of population

    etat01 = @SVector [0.1]
    etat02 = @SVector [1.0, 2.5]

    # modeling

    # Populations isolées
    mod_malthus = DynamiquePopulation.Modeling(etat01, [3, 2], (0.0, 20.0), 0.01, malthus)
    mod_logistic = DynamiquePopulation.Modeling(etat01, [1, 10], (0.0, 10.0), 0.01, logistic)
    mod_allee = DynamiquePopulation.Modeling(etat01, [1, 10, 2], (0.0, 3.0), 0.01, allee)
    
    # Populations en intéraction
    mod_lv = DynamiquePopulation.Modeling(etat02, [1.0, 1.0, 1.0, 1.0], (0.0, 30.0), 0.01, lv)
    mod_rma = DynamiquePopulation.Modeling(etat02, [1.0, 10.0, 1.0, 2.0, 2.0, 1.0], (0.0, 80.0), 0.01, rma)

    
    # test
    #simule(mod_lv, true, xlabel="Temps \$t\$", ylabel="Densités de population", label=["\$x(t)\$" "\$y(t)\$"], title="Simulation du modèle de Lotka Volterra ")
    #simule(mod_lv)

end

main()









###################################################source pour la création/structuration du module:##########################################################
###################################################### https://www.youtube.com/watch?v=WRI0ufhqPAU ##########################################################
#############################################################################################################################################################
