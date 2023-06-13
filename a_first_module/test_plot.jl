using Plots

function f(a::Int64, plot=true ; kwargs...)
    if plot
        Plots.plot(1:5,1:5; kwargs...)
    else
        a
    end
end

f(0, plot=false, title="test")