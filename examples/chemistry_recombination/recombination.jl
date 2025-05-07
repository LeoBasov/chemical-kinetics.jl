using ChemicalKinetics
using Plots
using LaTeXStrings

# simulation setup and execution
initialize!()

add_species!("data/O.json", mole_frac = 0.9)
add_species!("data/O2.json", mole_frac = 0.1)

add_reactions!("examples/chemistry_recombination/recombination.json")

set_T!(10000)
set_nrho!(1e23)
#set_relax_mode!("variable")

execute!(1e-4)

t, T = get_T(300)
t, T_O2 = get_Tvib(300, "O2")

fp_data = read_SPARTA_log("examples/chemistry_recombination/log.sparta")

t_fp = fp_data.dt * fp_data.data[1]["Step"]
T_fp = fp_data.data[1]["c_red_temp"]
Tvib_fp_O2 = fp_data.data[1]["c_red_tvib_O2"]
nrho_O_fp = fp_data.data[1]["c_red_nrho_O"]
nrho_O2_fp = fp_data.data[1]["c_red_nrho_O2"]

p = plot(t, T)
plot!(t, T_O2)
plot!(t_fp, T_fp)
plot!(t_fp, Tvib_fp_O2)

display(p)

t, nrho = get_nrho(300)

p = plot(t, nrho["O"], line = (2, :dashdot))
plot!(t, nrho["O2"], line = (2, :dashdot))
plot!(t_fp, nrho_O_fp)
plot!(t_fp, nrho_O2_fp)

display(p)