using ChemicalKinetics
using Plots
using Roots
using LaTeXStrings

initialize!()

add_species!("data/O2.json")

set_T!(10000)
set_Tvib!("O2", 7000)

set_nrho!(1e21)

set_relax_mode!("variable")

execute!(5e-4)

t, T = get_T(300)
t, Tvib_O2 = get_Tvib(300, "O2")

p = plot(t, T, line = 2, label="T - conti")
plot!(t, Tvib_O2, line = 2, label="Tvib O2 - conti")

display(p)

println("done")