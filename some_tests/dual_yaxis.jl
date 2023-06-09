using Plots

x = 1:0.1:10
y = x
z = x .^ 2


p1 = plot(x, y,
    label="P",
    legend=:topleft,
    xlabel="temps",
    ylabel="P",
    yguidefontrotation=-90,
    linestyle=:dashdotdot)
p1 = plot!(twinx(), x,z,
    c=:red,
    label="I",
    legend=:bottomright,
    ylabel="I",
    yguidefontrotation=-90,
    size=(400, 300),
    linestyle=:solid)


