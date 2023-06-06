using Plots

x = 1:0.1:10
y = x
z = x .^ 2


plot(x,y, label=false, ylabel="I")
plot!(twinx(), x,z,
    c=:red, label=false,
    ylabel="P",
    size=(400, 300))


