using ChemicalKinetics
using Plots

initialize!()

add_species!("data/NO.json", mole_frac = 0.2)
add_species!("data/N2.json", mole_frac = 0.2)
add_species!("data/N.json", mole_frac = 0.2)
add_species!("data/O2.json", mole_frac = 0.2)
add_species!("data/O.json", mole_frac = 0.2)

set_T!(10000)
set_nrho!(1e23)

set_Tvib!("NO", 10000)
set_Tvib!("N2", 10000)
set_Tvib!("O2", 10000)

read_reaction!("data/exchange.json")

execute!(4e-5)

t, T_NO = get_T(300, "NO")

p = plot(t, T_NO)

display(p)

t, X = get_molefrac(300)

display(plot(t, X))

println("done")