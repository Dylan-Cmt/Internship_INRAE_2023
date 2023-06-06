using Plots

x = 1:0.1:10
y = x
z = x .^ 2


p1 = plot(x, y,
    label=false,
    xlabel="rien",
    ylabel="P",
    yguidefontrotation=-90,
    linestyle=:dashdotdot)
p1 = plot!(twinx(), x,z,
    c=:red,
    label=false,
    ylabel="I",
    yguidefontrotation=-90,
    size=(400, 300),
    linestyle=:solid)


