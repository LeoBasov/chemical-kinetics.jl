using ChemicalKinetics

#add_species!("data/N2.json")
#add_species!("data/CH4.json")
#add_species!("data/CO2.json")

set_T!(300)
set_nrho!(1e22)

print_state()