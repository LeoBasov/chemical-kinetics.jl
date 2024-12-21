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

# simulation setup and execution
initialize!()

add_species!("data/NO.json", mole_frac = 0.2)
add_species!("data/N2.json", mole_frac = 0.2)
add_species!("data/N.json", mole_frac = 0.2)
add_species!("data/O2.json", mole_frac = 0.2)
add_species!("data/O.json", mole_frac = 0.2)

add_reactions!("data/exchange.json")

set_T!(10000)
set_nrho!(1e23)

execute!(4e-5)

# plotting
fp_data = read_log("examples/chemistry_exchange/chemistry_exchange.sparta")

t, T = get_T(300)
t, T_NO = get_Tvib(300, "NO")
t, T_N2 = get_Tvib(300, "N2")
t, T_O2 = get_Tvib(300, "O2")

p = plot(t, T, line = 2, label="conit")
#plot!(t, T_NO)
#plot!(t, T_N2)
#plot!(t, T_O2)

t_fp = fp_data[1] * 1e-8
T_fp = fp_data[2]
Tv_NO = fp_data[3]
Tv_N2 = fp_data[4]
Tv_O2 = fp_data[5]

plot!(t_fp, T_fp, line = (2, :dashdot), label="FP")
#plot!(t_fp, Tv_NO, line = (3, :dashdot))
#plot!(t_fp, Tv_N2, line = (3, :dashdot))
#plot!(t_fp, Tv_O2, line = (3, :dashdot))

xlabel!(L"t / s")
ylabel!(L"T / K")

display(p)

t, nrho = get_nrho(300)

nrho_NO = fp_data[6]
nrho_N2 = fp_data[7]
nrho_N = fp_data[8]
nrho_O2 = fp_data[9]
nrho_O = fp_data[10]

p = plot(t_fp, nrho_NO, line = (2, :dashdot), label="N2")
plot!(t_fp, nrho_N2, line = (2, :dashdot), label="O2")
plot!(t_fp, nrho_N, line = (2, :dashdot), label="O")
plot!(t_fp, nrho_O2, line = (2, :dashdot), label="NO")
plot!(t_fp, nrho_O, line = (2, :dashdot), label="N")

plot!(t, nrho[1], line = 2, label="NO")
plot!(t, nrho[2], line = 2, label="N2")
plot!(t, nrho[3], line = 2, label="N")
plot!(t, nrho[4], line = 2, label="O2")
plot!(t, nrho[5], line = 2, label="O")

xlabel!(L"t / s")
ylabel!(L"n_{\rho} / m^{-3}")

display(p)

println("done")