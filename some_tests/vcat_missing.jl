using Plots 

x=1:.1:10

y=vcat(x, missing, x.+10)
z = vcat(x, x .+ 10)

p1=plot(y.^2, y)
p2 = plot(z .^ 2, z)

# comparaison with/without the missing
plot(p1,p2,layout=(1,2))
