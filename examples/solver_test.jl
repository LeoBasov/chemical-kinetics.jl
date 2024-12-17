using DifferentialEquations
using Plots

function f(u, p, t)
    T = u[1] * 2.0 / (3.0*p.kb)
    eeq = p.kb*(3000) / (exp(3000/T) - 1)

    return [u[2] - u[1], u[1] - u[2]]
end

p = (kb=1.380649e-23, c=20)

tspan = (0.0, 10.0)
u0 = [5000, 1000]

prob = ODEProblem(f, u0, tspan, p);
sol = solve(prob);

plot(sol)

savefig("myplot.png")  

println("done")