using ChemicalKinetics
using Plots
using LaTeXStrings
using NCDatasets

kb = 1.380649e-23

# simulation setup and execution
initialize!()

add_species!("data/O.json", mole_frac = 0.9)
add_species!("data/O2.json", mole_frac = 0.1)

add_reactions!("examples/chemistry_recombination/recombination.json")

set_T!(10000)
set_nrho!(1e23)
set_relax_mode!("variable")

execute!(1e-3)

data_Georgii = NCDataset("examples/chemistry_recombination/iho__5June_O2rec_Arrhenius_10000K_VTVHScorrected_VTmolmolmixing.cdf")

t, T = get_T(300)
t, T_O2 = get_Tvib(300, "O2")

fp_data = read_SPARTA_log("examples/chemistry_recombination/log.sparta")

t_fp = fp_data.dt * fp_data.data[1]["Step"]
T_fp = fp_data.data[1]["c_red_temp"]
Tvib_fp_O2 = fp_data.data[1]["c_red_tvib_O2"]
nrho_O_fp = fp_data.data[1]["c_red_nrho_O"]
nrho_O2_fp = fp_data.data[1]["c_red_nrho_O2"]

ke_O = fp_data.data[1]["c_red_ke[1]"]
ke_O2 = fp_data.data[1]["c_red_ke[2]"]

erot_O = fp_data.data[1]["c_red_erot[1]"]
erot_O2 = fp_data.data[1]["c_red_erot[2]"]

evib_O = fp_data.data[1]["c_red_evib[1]"]
evib_O2 = fp_data.data[1]["c_red_evib[2]"]

p = plot(t, T, line = (2, :dashdot), label="conti T")
plot!(t, T_O2, line = (2, :dashdot), label="conti Tvib,O2")
plot!(t_fp, T_fp, label="FP T", line = 2)
plot!(t_fp, Tvib_fp_O2, label="FP Tvib,O2", line = 2)
plot!(t_fp, 20934.63283861*ones(length(t_fp)), label="FP Tvib,O2", line = 2)

plot!(data_Georgii["t"], data_Georgii["T"], label="Georgii T", line = (2, :dash))
plot!(data_Georgii["t"], data_Georgii["Tv_O2"], label="Georgii Tvib,O2", line = (2, :dash))

#ylims!(1000, 11000)
#xlims!(0, 1.5e-5)

display(p)

t, nrho = get_nrho(300)

#p = plot(t, nrho["O"], line = (2, :dashdot))
#plot!(t, nrho["O2"], line = (2, :dashdot))
p = plot(t_fp, nrho_O_fp, label="FP n_O", line = 2)
plot!(t_fp, nrho_O2_fp, label="FP n_O2", line = 2)
plot!(data_Georgii["t"], data_Georgii["n_O"], label="Georgii n_O", line = (2, :dash))
plot!(data_Georgii["t"], data_Georgii["n_O2"], label="Georgii n_O2", line = (2, :dash))

display(p)

ekin = nrho_O_fp.*ke_O + nrho_O2_fp.*ke_O2
erot = nrho_O_fp.*erot_O + nrho_O2_fp.*erot_O2
evib = nrho_O_fp.*evib_O + nrho_O2_fp.*evib_O2

p = plot(t_fp, ((ekin + erot + evib) / (ekin[begin] + erot[begin] + evib[begin]) - ones(length(ekin))))

xlabel!(L"t / s")
ylabel!(L"\frac{e_{tot} - e_{0}}{e_{0}}")

display(p)

println("done")