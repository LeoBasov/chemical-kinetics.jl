using Plots
using LaTeXStrings
using Revise
using Interpolations
using CSV
using DataFrames

function calc_E_O(T, nrho_O, kb)
    return 1.5 * kb * T .* nrho_O
end

function calc_E_Oe(T, nrho_O2, theta, kb)
    res = []

    for i in eachindex(T)
        push!(res, 1.5 * kb * T[i] * nrho_O2[i] + kb * T[i] * nrho_O2[i] + kb * theta / (exp(theta/T[i]) - 1.0) * nrho_O2[i])
    end

    return res
end

function cald_du_dt(t, u)
    du_dt = zeros(length(t))

    du_dt[begin] = (u[begin + 1] - u[begin]) / (t[begin + 1] - t[begin])
    du_dt[end] = (u[end] - u[end - 1]) / (t[end] - t[end - 1])

    for i in range(2, length(t) - 1)
        du_dt[i] = 0.5 * ((u[i + 1] - u[i]) / (t[i + 1] - t[i]) + (u[i] - u[i - 1]) / (t[i] - t[i - 1]))
    end

    return du_dt
end

const kb = 1.380649e-23
const mO = 2.65E-26
const mO2 = 5.31E-26
const theta = 2256.0
const N_A = 6.02214076e23

cantera_data = CSV.read("examples/dissociation_thermo_data/data/can.csv", DataFrame)
t = cantera_data.t/1000
T = cantera_data.T
nrho_O = cantera_data.yO .* cantera_data.rho / mO
nrho_O2 = cantera_data.yO2 .* cantera_data.rho / mO2

E_O = calc_E_O(T, nrho_O, kb)
E_O2 = calc_E_Oe(T, nrho_O2, theta, kb)

dnrho_O_dt = cald_du_dt(t, nrho_O)
dnrho_O2_dt = cald_du_dt(t, nrho_O2)
dE_O_dt = cald_du_dt(t, E_O)
dE_O2_dt = cald_du_dt(t, E_O2)

dE = dE_O_dt.*nrho_O + E_O.*dnrho_O_dt + dE_O2_dt.*nrho_O2 + E_O2.*dnrho_O2_dt

display(plot(T, dE./dnrho_O2_dt))

println("done")