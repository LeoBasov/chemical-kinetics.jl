using ChemicalKinetics
using Plots

initialize!()

add_species!("data/CH4.json", mole_frac = 1.0)

set_T!(10000)
set_nrho!(1e22)
set_Tvib!("CH4", [5000, 6000, 7000, 8000])

execute!(3e-5)

t, T = get_T(300)
t, Tvib = get_Tvib(300, "CH4")

p = plot(t, T)
plot!(t, Tvib)

display(p)

write2csv("CH4")

println("done")