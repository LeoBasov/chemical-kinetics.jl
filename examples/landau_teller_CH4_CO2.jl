using ChemicalKinetics
using Plots
using Roots

initialize!()

add_species!("data/CH4.json", mole_frac = 0.5)
add_species!("data/CO2.json", mole_frac = 0.5)

set_T!(10000)
set_nrho!(1e22)
set_Tvib!("CH4", 5000)
set_Tvib!("CO2", 4000)
set_Zvib!("CH4", [20, 30, 40, 50])
set_Zvib!("CO2", [20, 30, 40])

execute!(1e-5)

e0 = ChemicalKinetics.calc_etot(ChemicalKinetics._state)
f(T, p=ChemicalKinetics._state) = ChemicalKinetics.calc_etot(T, p) / e0 - 1.0

Teq = find_zero(f, 5000)

t, T = get_T(300)
t, T1 = get_Tvib(300, "CO2")
t, T2 = get_Tvib(300, "CH4")

Teq = ones(length(t)) * Teq

p = plot(t, T)
plot!(t, T1)
plot!(t, T2)
plot!(t, Teq, line = (3, :dashdot))

display(p)

write2csv("CH4_CO2")

println("done")