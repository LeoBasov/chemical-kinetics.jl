module ChemicalKinetics

export add_species!

include("Gas.jl")
include("Reader.jl")

state = State()

function add_species!(file_name)
    species = read_species(file_name)
    _add_species!(state, species)
end

end # module ChemicalKinetics
