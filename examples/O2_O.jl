using ChemicalKinetics
using Plots
using Roots

initialize!()

add_species!("data/O2.json", mole_frac = 0.5)
add_species!("data/O.json", mole_frac = 0.5)

set_T!(10000)
set_nrho!(1e22)
set_Tvib!("O2", 5000)

e0 = ChemicalKinetics.calc_etot(ChemicalKinetics._state)
f(T, p=ChemicalKinetics._state) = ChemicalKinetics.calc_etot(T, p) / e0 - 1.0

Teq = find_zero(f, 5000)

execute!(1e-2)

t, TO2 = get_T(300, "O2")

Teq = ones(length(t)) * Teq

plot(t, Teq, line = (3, :dashdot))
display(plot!(t, TO2))

t, X = get_molefrac(300)

display(plot(t, X))

println("done")