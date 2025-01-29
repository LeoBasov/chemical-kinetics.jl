using JSON

mutable struct SPARTAlog
    dt::Float64
    data::Array

    function SPARTAlog()
        new(0.0, [])
    end
end

function read_species(file_name)
    species = Species()
    
    open(file_name, "r") do file
        json = JSON.parse(file)

        species.name = json["name"]
        species.mass = json["mass"]
        species.dof_rot = json["rotation"]["dof"]
        species.dof_vib = json["vibration"]["dof"]
        species.Zrot = json["rotation"]["Z"]

        species.vhs.dref = json["vhs"]["dref"]
        species.vhs.Tref = json["vhs"]["Tref"]
        species.vhs.omega = json["vhs"]["omega"]
        species.vhs.alpha = json["vhs"]["alpha"]
        species.vhs.Zrotinf = json["vhs"]["Zrotinf"]
        species.vhs.Tstar = json["vhs"]["T*"]
        species.vhs.C1 = json["vhs"]["C1"]
        species.vhs.C2 = json["vhs"]["C2"]

        for mode in json["vibration"]["modes"]
            vibmode = Vibmode()

            vibmode.theta = mode["theta"]
            vibmode.degen = mode["degen"]
            vibmode.Z = mode["Z"]
            
            push!(species.vibmodes, vibmode)
        end
    end

    return species
end

function read_reactions(file_name)
    reactions = []
    
    open(file_name, "r") do file
        json = JSON.parse(file)
        
        for r in json["reactions"]
            reaction = Reaction()
            str = r["name"]
            splt = split(replace(str,"+"=>""), " ")
            pre = true

            for elem in splt
                if elem == "-->"
                    pre = false
                elseif pre ==true && !(elem == "")
                    if elem in keys(reaction.stochio_coeff)
                        reaction.stochio_coeff[elem] -= 1
                    else
                        reaction.stochio_coeff[elem] = -1
                        push!(reaction.reactants, elem)
                    end
                elseif pre ==false && !(elem == "")
                    if elem in keys(reaction.stochio_coeff)
                        reaction.stochio_coeff[elem] += 1
                    else
                        reaction.stochio_coeff[elem] = 1
                    end
                end
            end

            reaction.A = r["Arrhenius"]["A"]
            reaction.B = r["Arrhenius"]["B"]
            reaction.Ea = r["Arrhenius"]["Ea"]
            reaction.DeltaE = r["Arrhenius"]["DeltaE"]

            push!(reactions, reaction)
        end
    end

    return reactions
end

function _setup_data(splt)
    data = Dict()
    names = []
    
    for elem in splt
        if elem != ""
            push!(names, elem)
            data[elem] = []
        end
    end
    
    return names, data
end

function _read_data(file_name)
    found = false
    read = false
    data_loc = Dict()
    data = []
    names = []

    open(file_name, "r") do file
        for line in readlines(file)
            splt = split(line, " ")
            len = length(splt)
            
            if len > 1 && splt[1] == "Step"
                found = true
                data_loc = Dict()
            end

            if found == true && read == false
                names, data_loc = _setup_data(splt)
                read = true
                continue
            end
            
            if found && read && length(splt) != length(names)
                found = false
                read = false
                push!(data, data_loc)
            end
            
            if found == true && read == true
                try
                    for i in eachindex(names)
                        push!(data_loc[names[i]], parse(Float64, splt[i]))
                    end
                catch err
                    break
                end
            end
        end
    end

    return data
end

function _read_timestep(file_path)
    dt = -1.0

    open(file_path, "r") do file
        for line in readlines(file)
            splt = split(line, " ")
            len = length(splt)
    
            if len > 1 && splt[1] == "timestep"
                dt = parse(Float64, splt[end])
                break
            end
        end
    end

    return dt
end

function read_SPARTA_log(file_path)
    log = SPARTAlog()

    log.dt = _read_timestep(file_path)
    log.data = _read_data(file_path)

    return log
end