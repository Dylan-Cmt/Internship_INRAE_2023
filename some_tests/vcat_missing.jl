using Plots 

x=1:.1:10

y=vcat(x, missing, x.+10)

plot(y.^2, y)

