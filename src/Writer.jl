function write2csv(file_prefix, N = 1000)
    println("writing " * file_prefix * ".csv")

    open(file_prefix * ".csv", "w") do file
        str = "t,T"
        str_Tvib = ""
        str_nrho = ""
        t, T = get_T(N)
        t, nrho = get_nrho(N)
        Tvib = Dict()

        for species in _state.species
            str_nrho *= ",nrho_" * species.first
            t, Tvib_s = get_Tvib(N, species.first)
            Tvib[species.first] = Tvib_s

            for v in 1:length(species.second.vibmodes)
                str_Tvib *= ",Tvib_" * species.first * "_" * string(v)
            end
        end

        write(file, str * str_Tvib * str_nrho * "\n")

        for i in 1:N
            str = string(t[i]) * "," * string(T[i])

            for species in _state.species
                for v in 1:length(species.second.vibmodes)
                    str *= "," * string(Tvib[species.first][v][i])
                end
            end

            for k in 1:length(_state.species)
                str *= "," * string(nrho[k][i])
            end

            write(file, str * "\n")
        end
    end
end