using ChemicalKinetics
using Plots
using Roots
using LaTeXStrings

function read_log(file_name)
    data = [[], [], [], [], [], [], [], [], [], []]

    open(file_name, "r") do file
        for line in readlines(file)
            splt = split(line, " ")
            filter!(e->e≠"",splt)
            push!(data[1], parse(Float64, splt[1]))
            push!(data[2], parse(Float64, splt[3]))
            push!(data[3], parse(Float64, splt[5]))
            push!(data[4], parse(Float64, splt[6]))
            push!(data[5], parse(Float64, splt[7]))

            push!(data[6], parse(Float64, splt[8]))
            push!(data[7], parse(Float64, splt[9]))
            push!(data[8], parse(Float64, splt[10]))
            push!(data[9], parse(Float64, splt[11]))
            push!(data[10], parse(Float64, splt[12]))
        end
    end

    return data
end

initialize!()

add_species!("data/NO.json", mole_frac = 0.2)
add_species!("data/N2.json", mole_frac = 0.2)
add_species!("data/N.json", mole_frac = 0.2)
add_species!("data/O2.json", mole_frac = 0.2)
add_species!("data/O.json", mole_frac = 0.2)

set_T!(10000)
set_nrho!(1e23)

set_Tvib!("NO", 5000)
set_Tvib!("N2", 5000)
set_Tvib!("O2", 5000)

execute!(1e-3)

e0 = ChemicalKinetics.calc_etot(ChemicalKinetics._state)
f(T, p=ChemicalKinetics._state) = ChemicalKinetics.calc_etot(T, p) / e0 - 1.0

Teq = find_zero(f, 5000)

fp_data = read_log("examples/thermal_relaxation/therm_relax_log.sparta")

t, T = get_T(300)
t, T_NO = get_Tvib(300, "NO")
t, T_N2 = get_Tvib(300, "N2")
t, T_O2 = get_Tvib(300, "O2")

p = plot(t, T, line = 2, label="T - conti")
plot!(t, T_NO, line = 2, label="Tvib NO - conti")
plot!(t, T_N2, line = 2, label="Tvib N2 - conti")
plot!(t, T_O2, line = 2, label="Tvib O2 - conti")

Teq = ones(length(t)) * Teq

plot!(t, Teq, line = (2, :dash), label="Teq")

t_fp = fp_data[1] * 1e-7
T_fp = fp_data[2]
Tv_NO = fp_data[3]
Tv_N2 = fp_data[4]
Tv_O2 = fp_data[5]

plot!(t_fp, T_fp, line = (2, :dashdot), label="T - FP")
plot!(t_fp, Tv_NO, line = (2, :dashdot), label="Tvib NO - FP")
plot!(t_fp, Tv_N2, line = (2, :dashdot), label="Tvib N2 - FP")
plot!(t_fp, Tv_O2, line = (2, :dashdot), label="Tvib O2 - FP")

xlabel!(L"t / s")
ylabel!(L"T / K")

display(p)

println("done")