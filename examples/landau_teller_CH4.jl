using ChemicalKinetics
using Plots

add_species!("data/CH4.json", mole_frac = 1.0)

set_T!(10000)
set_nrho!(1e22)
set_Tvib!("CH4", 5000)

solve!(1e-5)

#sol = get_energy(10)



t, T = get_T(300)

"""sol = ChemicalKinetics._solution

plot(sol.t, sol[1, :])
plot!(sol.t, sol[2, :])
plot!(sol.t, sol[3, :])
plot!(sol.t, sol[4, :])"""

display(plot(t, T))

println("done")