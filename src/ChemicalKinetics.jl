module ChemicalKinetics

export set_nrho!
export set_T!
export set_frac!
export add_species!
export print_state
export clear!

include("Gas.jl")
include("Reader.jl")

state = State()

function clear!()
    global state = State()
end

function set_nrho!(nrho)
    state.nrho = nrho
end

function set_T!(T)
    state.T = T
end

function set_frac!(species, mole_frac)
    state.mole_fractions[species] = mole_frac
end

function add_species!(file_name)
    species = read_species(file_name)
    _add_species!(state, species)
end

function print_state()
    println("T: ", state.T)
    println("nrho: ", state.nrho)
    println("mole fractions:")

    for key in keys(state.mole_fractions)
        println(key, ": ", state.mole_fractions[key])
    end
end

end # module ChemicalKinetics
