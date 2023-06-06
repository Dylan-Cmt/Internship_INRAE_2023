using Plots

x = 1:0.1:10
y = x
z = x .^ 2


p1 = plot(x, y,
    label=false,
    ylabel="P",
    linestyle=:dashdotdot)
p1 = plot!(twinx(), x,z,
    c=:red,
    label=false,
    ylabel="I",
    size=(400, 300),
    linestyle=:solid)


