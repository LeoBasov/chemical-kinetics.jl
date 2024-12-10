using ChemicalKinetics

add_species!("data/CH4.json", mole_frac = 1.0)

set_T!(300)
set_nrho!(1e22)
set_Tvib!("CH4", [100, 200, 300, 400])

print_state()