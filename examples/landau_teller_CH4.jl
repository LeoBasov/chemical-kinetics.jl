using ChemicalKinetics
using Plots
using NCDatasets

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

#=
write2netCDF("CH4")

ds = NCDataset("CH4.nc")

t = ds["time"]
T = ds["temperature"]
nrho_CH4 = ds["nrho_CH4"]
Tvib_CH4_1 = ds["Tvib_CH4_1"]
Tvib_CH4_2 = ds["Tvib_CH4_2"]
Tvib_CH4_3 = ds["Tvib_CH4_3"]
Tvib_CH4_4 = ds["Tvib_CH4_4"]

p = plot(t, T)
p = plot!(t, Tvib_CH4_1)
p = plot!(t, Tvib_CH4_2)
p = plot!(t, Tvib_CH4_3)
p = plot!(t, Tvib_CH4_4)
display(p)

println("done")
=#