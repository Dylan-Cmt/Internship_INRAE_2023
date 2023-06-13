using Plots

x = 1:1:365*2
y = cos.(x/50)

vspan([0, 365 / 2, 365, 3 * 365 / 2, 2 * 365],
    color=:lightgray,
    label="growing season")

# vspan!([365 / 2, 365, 2*365/3, 2*365],
#     color=:white,
#     label="winter season")
    
plot!(x,y,
    label=false,
    c=:black,
    xlabel="x",
    ylabel="y")