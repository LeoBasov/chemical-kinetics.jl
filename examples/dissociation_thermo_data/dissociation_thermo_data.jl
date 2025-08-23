using ChemicalKinetics

# simulation setup and execution
initialize!()

add_species!("data/O.json", mole_frac = 0.2)
add_species!("data/O2.json", mole_frac = 0.8)

add_reactions!("examples/dissociation_thermo_data/reactions.json")

set_T!(10000)
set_nrho!(1e23)
set_relax_mode!("variable")

execute!(3e-4)

t, T = get_T(300)
t, Tvib = get_Tvib(300, "O2")
t, nrho = get_nrho(300)

p = plot(t, T, label=L"T")
plot!(t, Tvib, label=L"T_{\mathrm{vib}}")

xlabel!(L"t\,/\,\mathrm{s}")
ylabel!(L"T\,/\,\mathrm{K}")

display(p)

p = plot(t, nrho["O2"], label=L"\mathrm{O}_2")
plot!(t, nrho["O"], label=L"\mathrm{O}")

xlabel!(L"t\,/\,\mathrm{s}")
ylabel!(L"n\,/\,\mathrm{m}^{-3}")

display(p)

println("done")