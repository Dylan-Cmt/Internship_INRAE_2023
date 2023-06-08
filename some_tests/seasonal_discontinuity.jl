using Plots 

x1=1:.1:5
x2=10:.1:15
x3 = 20:.1:30

y1 = x1
y2 = sqrt.(x2) 
y3 = x3 .^(1/3)

plot([x1, x2, x3],[y1, y2, y3])