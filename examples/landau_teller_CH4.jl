using ChemicalKinetics
using Plots

add_species!("data/CH4.json", mole_frac = 1.0)

set_T!(10000)
set_nrho!(1e22)
set_Tvib!("CH4", 5000)

sol = execute(1.0)

plot(sol);
plot!()