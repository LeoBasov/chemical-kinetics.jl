using ChemicalKinetics
using Plots
using LaTeXStrings
using Revise

# simulation setup and execution
initialize!()

add_species!("examples/dissociation_thermo_data/data/O.json", mole_frac = 0.2)
add_species!("examples/dissociation_thermo_data/data/O2.json", mole_frac = 0.8)

add_reactions!("examples/dissociation_thermo_data/data/reactions.json")

set_T!(10000)
set_nrho!(1e23)
#set_relax_mode!("variable")
set_reaction_enthalpy_mode!("constant")

t_max = 1e-4
execute!(t_max)

t, T = get_T(300)
t, Tvib = get_Tvib(300, "O2")
t, nrho = get_nrho(300)

cantera_data = read_csv("examples/dissociation_thermo_data/data/can.csv")

p = plot(t, T, label=L"T - \mathrm{conti}")
plot!(t, Tvib, label=L"T_{\mathrm{vib}} - \mathrm{conti}")
plot!(cantera_data.t/1000, cantera_data.T, label=L"T - \mathrm{Cantera}")

xlabel!(L"t\,/\,\mathrm{s}")
ylabel!(L"T\,/\,\mathrm{K}", xlim=(0, t_max))

display(p)

p = plot(t, nrho["O2"], label=L"\mathrm{O}_2 - \mathrm{conti}")
plot!(t, nrho["O"], label=L"\mathrm{O} - \mathrm{conti}")
plot!(cantera_data.t/1000, cantera_data.yO2.*cantera_data.rho/5.31E-26, label=L"\mathrm{O}_2 - \mathrm{Cantera}", line = (2, :dash))
plot!(cantera_data.t/1000, cantera_data.yO.*cantera_data.rho/2.65E-26, label=L"\mathrm{O} - \mathrm{Cantera}", xlim=(0, t_max), line = (2, :dash))


xlabel!(L"t\,/\,\mathrm{s}")
ylabel!(L"n\,/\,\mathrm{m}^{-3}")

display(p)

println("done")