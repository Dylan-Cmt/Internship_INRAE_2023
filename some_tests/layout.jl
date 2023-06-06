using Plots

x = 1:0.1:10
y = x
z = x .^ 2

p1 = plot(x,y)
p2 = plot(x,z)
plot(p1,p2,layout=(2,1))