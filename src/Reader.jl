using JSON

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
        reaction = Reaction()
    end

    return reactions
end