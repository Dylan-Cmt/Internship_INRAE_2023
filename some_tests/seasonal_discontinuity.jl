using Plots 

x1 = 1:.1:5
x2 = 5:.1:15
x3 = 15:.1:30
x4 = 40:.1:45

y1 = x1
y2 = sqrt.(x2) 
y3 = x3 .^(1/3)
y4 = x4 .^(1/4)


t = [[]]                    # initialize a vector of vectors
push!(t, x1)                # add each vector of time
push!(t, x2)    
push!(t, x3)
push!(t,x4)
popfirst!(t)                # pop the empty first

x = [[]]                    # again
push!(x,y1)
push!(x, y2)
push!(x, y3)
push!(x,y4)
popfirst!(x)

plot(t,x)