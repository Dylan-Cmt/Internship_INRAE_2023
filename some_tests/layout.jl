using Plots

# This can be useful if we want to plot both model in the same window

# mod 1
x = 1:0.1:10
y = x
z = x .^ 2

# mod 2
y1 = x .- 1
z1 = x .^ 2 .- 1

# create top plot
p1 = plot(x,y)
# add bot model
p1 = plot!(x,y1)
# create downside plot
p2 = plot(x,z)
# add other model
p2 = plot!(x,z1)

# plot the superposed models in the same window
plot(p1,p2,layout=(2,1))