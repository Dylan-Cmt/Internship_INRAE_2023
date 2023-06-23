using Plots

x = 1:.1:10
y = x.^2
y_err = 1

plot(x, y, ribbon=(20, 10), label="Data with Uncertainty")
