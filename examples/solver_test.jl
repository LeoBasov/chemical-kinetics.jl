using DifferentialEquations
using Plots

function f(u, p, t)
    return u
end

tspan = (0.0, 1.0)
u0 = 1.0

prob = ODEProblem(f, u0, tspan);
sol = solve(prob);

plot(sol)

savefig("myplot.png")  

println("done")