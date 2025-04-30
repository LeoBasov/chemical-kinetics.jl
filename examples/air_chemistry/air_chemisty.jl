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

add_reactions!("examples/air_chemistry/air_chemistry.json")

set_T!(10000)
set_nrho!(1e23)
set_relax_mode!("variable")

execute!(1e-5)

# plotting
fp_data = read_SPARTA_log("examples/air_chemistry/log.sparta")

t, T = get_T(300)
t, T_NO = get_Tvib(300, "NO")
t, T_N2 = get_Tvib(300, "N2")
t, T_O2 = get_Tvib(300, "O2")

p = plot(t, T, line = 2, label="conit")
#plot!(t, T_NO)
#plot!(t, T_N2)
#plot!(t, T_O2)

t_fp = fp_data.dt * fp_data.data[1]["Step"]
T_fp = fp_data.data[1]["c_red_temp"]
#Tv_NO = fp_data[3]
#Tv_N2 = fp_data[4]
#Tv_O2 = fp_data[5]

plot!(t_fp, T_fp, line = (2, :dashdot), label="FP")
#plot!(t_fp, Tv_NO, line = (3, :dashdot))
#plot!(t_fp, Tv_N2, line = (3, :dashdot))
#plot!(t_fp, Tv_O2, line = (3, :dashdot))

xlabel!(L"t / s")
ylabel!(L"T / K")

display(p)

t, nrho = get_nrho(300)

nrho_NO = fp_data.data[1]["c_red_nrho_NO"]
nrho_N2 = fp_data.data[1]["c_red_nrho_N2"]
nrho_N = fp_data.data[1]["c_red_nrho_N"]
nrho_O2 = fp_data.data[1]["c_red_nrho_O2"]
nrho_O = fp_data.data[1]["c_red_nrho_O"]

p = plot(t_fp, nrho_NO, line = (2, :dashdot), label="N2")
plot!(t_fp, nrho_N2, line = (2, :dashdot), label="O2")
plot!(t_fp, nrho_N, line = (2, :dashdot), label="O")
plot!(t_fp, nrho_O2, line = (2, :dashdot), label="NO")
plot!(t_fp, nrho_O, line = (2, :dashdot), label="N")

plot!(t, nrho["NO"], line = 2, label="NO")
plot!(t, nrho["N2"], line = 2, label="N2")
plot!(t, nrho["N"], line = 2, label="N")
plot!(t, nrho["O2"], line = 2, label="O2")
plot!(t, nrho["O"], line = 2, label="O")

xlabel!(L"t / s")
ylabel!(L"n_{\rho} / m^{-3}")

xlims!(0, 2e-6)

display(p)

println("done")