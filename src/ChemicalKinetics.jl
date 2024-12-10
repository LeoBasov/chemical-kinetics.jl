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

_state::State = State()
_verbose::Bool = true

function initialize!(;verbose::Bool = true)
    global _state = State()
    global _verbose = verbose
end

function  set_verbosity!(verbose)
    global _verbose = verbose
end

function set_nrho!(nrho)
    _state.nrho = nrho
end

function set_T!(T)
    _state.T = T
end

function set_molefrac!(species_name, mole_frac)
    _state.mole_fractions[species_name] = mole_frac
    if _verbose == true
        println("set mole fraction of [" * species_name * "] to: " * string(mole_frac))
    end
end

function add_species!(file_name; mole_frac = 0.0)
    species = read_species(file_name)
    _add_species!(_state, species)
    if _verbose == true
        println("added species [" * species.name * "]")
    end
    set_molefrac!(species.name, mole_frac)
end

function print_state()
    println("T: ", _state.T)
    println("nrho: ", _state.nrho)
    println("mole fractions:")

    for key in keys(_state.mole_fractions)
        println(key, ": ", _state.mole_fractions[key])
    end
end

end # module ChemicalKinetics
