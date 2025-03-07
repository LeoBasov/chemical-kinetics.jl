using ChemicalKinetics
using Plots
using Roots
using LaTeXStrings

initialize!()

add_species!("data/O2.json", mole_frac = 0.5)
add_species!("data/N2.json", mole_frac = 0.5)

set_T!(10000)
set_Tvib!("O2", 7000)
set_Tvib!("N2", 5000)

set_nrho!(2e21)

set_relax_mode!("variable")

execute!(1e-3)

# RESULTS
sparta_log = read_SPARTA_log("examples/O2_N2_thermal_relaxation/log.sparta")
t_sp = sparta_log.dt * sparta_log.data[1]["Step"]
T_sp = sparta_log.data[1]["c_temp_red"]
Tvib_O2_sp = sparta_log.data[1]["c_tvib_O2_red"]
Tvib_N2_sp = sparta_log.data[1]["c_tvib_N2_red"]

t, T = get_T(300)
t, Tvib_O2 = get_Tvib(300, "O2")
t, Tvib_N2 = get_Tvib(300, "N2")

p = plot(t, T, line = 2, label="T - conti")
plot!(t, Tvib_O2, line = 2, label="Tvib O2 - conti")
plot!(t, Tvib_N2, line = 2, label="Tvib N2 - conti")

plot!(t_sp, T_sp, line = (2, :dashdot), label="T - FP")
plot!(t_sp, Tvib_O2_sp, line = (2, :dashdot), label="Tvib O2 - FP")
plot!(t_sp, Tvib_N2_sp, line = (2, :dashdot), label="Tvib N2 - FP")

display(p)

println("done")