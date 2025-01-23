using NCDatasets

function write2netCDF(file_prefix, N = 1000)
    println("writing " * file_prefix * ".nc")

    NCDataset(file_prefix * ".nc","c") do ds
        defDim(ds, "length", N)
        defDim(ds, "Nspecies", length(_state.species))

        t, T = get_T(N)
        t, nrho = get_nrho(N)
        k = 1

        time_var = defVar(ds, "time", Float64, ("length",))
        temp_var = defVar(ds, "temperature", Float64, ("length",))

        time_var[:] = t
        temp_var[:] = T

        for species in _state.species
            nrho_var = defVar(ds, "nrho_" * species.first, Float64, ("length", ))
            nrho_var[:] = nrho[k]
            k += 1

            t, Tvib_s = get_Tvib(N, species.first)
            j = 1

            for v in 1:length(species.second.vibmodes)
                str_Tvib = "Tvib_" * species.first * "_" * string(v)
                Tvib_var = defVar(ds, str_Tvib, Float64, ("length", ))
                Tvib_var[:] = Tvib_s[j]
                j += 1
            end
        end
        
    end
end

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