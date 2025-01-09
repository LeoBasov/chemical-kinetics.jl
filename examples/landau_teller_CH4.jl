using ChemicalKinetics
using Plots

initialize!()

add_species!("data/CH4.json", mole_frac = 1.0)

set_T!(10000)
set_nrho!(1e22)
set_Tvib!("CH4", 5000)
set_Zvib!("CH4", [20, 30, 40, 50])

execute!(1e-5)

t, T = get_T(300)
t, Tvib = get_Tvib(300, "CH4")

p = plot(t, T)
plot!(t, Tvib)

display(p)

write2csv("CH4")

println("done")