using Plots

τ = 120
Τ = 365
n = 5

x = 1:1:Τ*n
y = cos.(x / 50)

v = [[0, τ]]
vspan(v[1] ./ Τ,
    color=:lightgray,
    label="growing season")

for i in 1:n-1
    u = [i * Τ, i * Τ + τ]
    push!(v, u)
end



vspan!(v[2:end] ./ Τ,
    color=:lightgray,
    label=false)

annee = x ./ Τ
plot!(annee, y,
    label=false,
    c=:black,
    xlabel="x",
    ylabel="y")