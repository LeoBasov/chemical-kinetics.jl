using ChemicalKinetics
using Plots
using Roots
using LaTeXStrings

function particle(nuN, nuN2, nN0, nN20, k, V, w, dt, tmax)
    time = [0.0]
    nN = [Int(nN0 * V / w)]
    nN2 = [Int(nN20 * V / w)]
    A = []

    println("N - N: ", nN[1])
    println("N - N2: ", nN2[1])

    for i in 1:Int(round(tmax/dt))
        Nmin = min(nN[end], nN2[end])
        Nmax = max(nN[end], nN2[end])
        a = Nmax*w*k / V
        R = rand(Nmin)
        p = 1.0 - exp(-a*dt)

        push!(time, time[end] + dt)
        push!(nN, nN[end])
        push!(nN2, nN2[end])
        push!(A, a)

        for j in 1:Nmin
            if R[j] < p
                nN[end] += nuN
                nN2[end] += nuN2
            end
        end
    end

    return time, nN*w/V, nN2*w/V, A
end

function analytical_nN2(nN, nN2, k, t)
    ntot = nN + nN2
    Q = ntot + nN2
    num = nN2 * Q
    denum = 2 * nN2 - (nN2 - ntot) * exp(Q*k*t)

    return num / denum
end

function get_sol(nN, nN2, k, tmax, N = 100)
    dt = tmax/N
    t = [0.0]
    solN = [nN]
    solN2 = [nN2]

    for i in 2:N
        tt = t[i - 1] + dt
        nnN2 = analytical_nN2(nN, nN2, k, tt)
        nnN = (nN2 - nnN2)*2 + nN
        push!(t, tt)
        push!(solN, nnN)
        push!(solN2, nnN2)
    end

    return t, solN, solN2
end

# simulation setup and execution
initialize!()

ntot = 1e21
XN = 0.02
XN2 = 0.98
tmax = 5e-7
k = 1e-14
V = 0.1
dt = 5e-10
w = 1e14

add_species!("data/N2.json", mole_frac = XN2)
add_species!("data/N.json", mole_frac = XN)

add_reactions!("examples/dissociation_analytical/reaction.json")

set_T!(10000)
set_nrho!(ntot)

execute!(tmax)

tana, solN, solN2 = get_sol(ntot*XN, ntot*XN2, 1e-14, tmax)
time, nN, nN2, A = particle(2, -1, XN*ntot, XN2*ntot, k, V, w, dt, tmax)

t, T = get_T(300)
t, Tvib = get_Tvib(300, "N2")
t, nrho = get_nrho(300)

p = plot(tana, solN2, line = (2, :dashdot), label="N2 - analytic")
plot!(tana, solN, line = (2, :dashdot), label="N - analytic")
plot!(t, nrho["N"], line=2, label="N")
plot!(t, nrho["N2"], line=2, label="N2")
plot!(time, nN, line = (2, :dash), label="N - particle")
plot!(time, nN2, line = (2, :dash), label="N2 - particle")
display(p)

println("dt: ", minimum([1e-2/a for a in A]))