module EDOs

export malthus, logistic, allee, lv, rma


# Imports nécessaires
using StaticArrays                                  # pour l'utilisation de @SVector 


# Modèle de Malthus
function malthus(u, p, t)
    n, m = p                                        # unpack the vectors into scalar
    x = u[1]
    dx = (n - m) * x                                # dot x
    @SVector [dx]                                   # return a new vector
end

# Modèle Logistique
function logistic(u, p, t)
    r, K = p                                        # unpack the vectors into scalar
    x = u[1]
    dx = r * x * (1 - x / K)                        # dot x
    @SVector [dx]                                   # return a new vector
end

# Modèle de l'effet Allee
function allee(u, p, t)
    r, K, ϵ = p                                     # unpack the vectors into scalar
    x = u[1]
    dx = r * x * (1 - x / K) * (x / ϵ - 1)
    @SVector [dx]                                   # return a new vector
end

# Modèle de Lotka-Volterra
function lv(u, params, t)
    r, c, b, m = params                              # unpack the vectors into scalar
    x = u[1]
    y = u[2]
    dx = r * x - c * x * y                           # dot x
    dy = b * x * y - m * y                           # dot y
    @SVector [dx, dy]                                # return a new vector
end

# Modèle de Rosenzweig-MacArthur
function rma(u, params, t)
    r, K, c, h, b, m = params                        # unpack the vectors into scalar
    x = u[1]
    y = u[2]
    dx = r * x * (1 - x / K) - c * x / (h + x) * y   # dot x
    dy = b * x / (h + x) * y - m * y                 # dot y
    @SVector [dx, dy]                                # return a new vector
end

end # end of the module