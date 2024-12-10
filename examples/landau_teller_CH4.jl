using ChemicalKinetics
using Plots

add_species!("data/N2.json", mole_frac = 1.0)

set_T!(10000)
set_nrho!(1e22)
set_Tvib!("N2", 5000)

sol = execute(2e-2)

plot(sol)
plot!()