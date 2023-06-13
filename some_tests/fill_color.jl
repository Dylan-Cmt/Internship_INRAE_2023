using Plots

x = collect(range(0, 2, length=100))
y1 = exp.(x)
y2 = exp.(1.3 .* x)


plot(x, y1, fillrange=y2, fillalpha=0.2, color=:lightgray, label="winter", legend=:topleft)