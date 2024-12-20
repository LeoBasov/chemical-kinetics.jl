using ChemicalKinetics
using Plots
using Roots

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

fp_data = read_log("examples/chemistry_exchange/chemistry_exchange.sparta")

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
t, T_N2 = get_T(300, "N2")
t, T_O2 = get_T(300, "O2")

p = plot(t, T_NO)
p = plot!(t, T_N2)
p = plot!(t, T_O2)

t_fp = fp_data[1] * 1e-8
T_fp = fp_data[2]
Tv_NO = fp_data[3]
Tv_N2 = fp_data[4]
Tv_O2 = fp_data[5]

p = plot!(t_fp, T_fp, line = (3, :dashdot))
p = plot!(t_fp, Tv_NO, line = (3, :dashdot))
p = plot!(t_fp, Tv_N2, line = (3, :dashdot))
p = plot!(t_fp, Tv_O2, line = (3, :dashdot))

display(p)

t, X = get_molefrac(300)

nrho_NO = fp_data[6]
nrho_N2 = fp_data[7]
nrho_N = fp_data[8]
nrho_O2 = fp_data[9]
nrho_O = fp_data[10]

plot(t_fp, nrho_NO / 1e23, line = (3, :dashdot))
plot!(t_fp, nrho_N2 / 1e23, line = (3, :dashdot))
plot!(t_fp, nrho_N / 1e23, line = (3, :dashdot))
plot!(t_fp, nrho_O2 / 1e23, line = (3, :dashdot))
plot!(t_fp, nrho_O / 1e23, line = (3, :dashdot))

display(plot!(t, X))

println("done")