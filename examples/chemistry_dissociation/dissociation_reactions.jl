using ChemicalKinetics
using Plots
using Roots
using LaTeXStrings

# simulation setup and execution
initialize!()

add_species!("data/NO.json", mole_frac = 0.2)
add_species!("data/N2.json", mole_frac = 0.2)
add_species!("data/N.json", mole_frac = 0.2)
add_species!("data/O2.json", mole_frac = 0.2)
add_species!("data/O.json", mole_frac = 0.2)

add_reactions!("data/dissociation.json")

set_T!(10000)
set_nrho!(1e23)

execute!(4e-5)

# plotting
log = read_SPARTA_log("examples/chemistry_dissociation/dissociation.sparta")

t, T = get_T(300)
t, T_NO = get_Tvib(300, "NO")
t, T_N2 = get_Tvib(300, "N2")
t, T_O2 = get_Tvib(300, "O2")

p = plot(t, T, line = 2, label="T - conit", linecolor=:black)
plot!(t, T_NO, line = 2, label="Tv NO - conit", linecolor=:blue)
plot!(t, T_N2, line = 2, label="Tv N2 - conit", linecolor=:green)
plot!(t, T_O2, line = 2, label="Tv O2 - conit", linecolor=:red)

t_fp = log.dt * log.data[1]["Step"]
T_fp = log.data[1]["c_red_temp"]
Tv_NO = log.data[1]["c_red_tvib_NO"]
Tv_N2 = log.data[1]["c_red_tvib_N2"]
Tv_O2 = log.data[1]["c_red_tvib_O2"]

plot!(t_fp, T_fp, line = (2, :dashdot), label="T - FP", linecolor=:black)
plot!(t_fp, Tv_NO, line = (2, :dashdot), label="Tv NO - FP", linecolor=:blue)
plot!(t_fp, Tv_N2, line = (2, :dashdot), label="Tv N2 - FP", linecolor=:green)
plot!(t_fp, Tv_O2, line = (2, :dashdot), label="Tv O2 - FP", linecolor=:red)

xlabel!(L"t / s")
ylabel!(L"T / K")

display(p)

t, nrho = get_nrho(300)

nrho_NO = log.data[1]["c_red_nrho_NO"]
nrho_N2 = log.data[1]["c_red_nrho_N2"]
nrho_N = log.data[1]["c_red_nrho_N"]
nrho_O2 = log.data[1]["c_red_nrho_O2"]
nrho_O = log.data[1]["c_red_nrho_O"]

p = plot(t_fp, nrho_NO, line = (2, :dashdot), label="N2")
plot!(t_fp, nrho_N2, line = (2, :dashdot), label="O2")
plot!(t_fp, nrho_N, line = (2, :dashdot), label="O")
plot!(t_fp, nrho_O2, line = (2, :dashdot), label="NO")
plot!(t_fp, nrho_O, line = (2, :dashdot), label="N")

plot!(t, nrho[1], line = 2, label="NO")
plot!(t, nrho[2], line = 2, label="N2")
plot!(t, nrho[3], line = 2, label="N")
plot!(t, nrho[4], line = 2, label="O2")
plot!(t, nrho[5], line = 2, label="O")

xlabel!(L"t / s")
ylabel!(L"n_{\rho} / m^{-3}")

display(p)

println("done")