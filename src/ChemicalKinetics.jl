module ChemicalKinetics

export set_nrho!
export set_T!
export set_molefrac!
export add_species!
export print_state
export initialize!

include("Gas.jl")
include("Reader.jl")
include("Solver.jl")
include("Constants.jl")

state::State = State()
_verbose::Bool = true

function initialize!(;verbose::Bool = true)
    global state = State()
    global _verbose = verbose
end

function  set_verbosity!(verbose)
    global _verbose = verbose
end

function set_nrho!(nrho)
    state.nrho = nrho
end

function set_T!(T)
    state.T = T
end

function set_molefrac!(species_name, mole_frac)
    state.mole_fractions[species_name] = mole_frac
    if _verbose == true
        println("set mole fraction of [" * species_name * "] to: " * string(mole_frac))
    end
end

function add_species!(file_name; mole_frac = 0.0)
    species = read_species(file_name)
    _add_species!(state, species)
    if _verbose == true
        println("added species [" * species.name * "]")
    end
    set_molefrac!(species.name, mole_frac)
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
