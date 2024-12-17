module ChemicalKinetics

export set_nrho!
export set_T!
export set_Tvib!
export set_molefrac!
export add_species!
export print_state
export initialize!
export solve!
export get_energy
export get_T
export get_molefrac

include("Gas.jl")
include("Reader.jl")
include("Solver.jl")
include("Constants.jl")

_state::State = State()
_verbose::Bool = true
_solution = nothing
_tmax::Number = 0.0

function get_molefrac(N)
    t = []
    X = []

    for i in 1:length(_state.species)
        push!(X, [])
    end

    for tt in range(0, _tmax, N)
        push!(t, tt * t_tilde)
        i = 1

        for species in _state.species
            molde_fraction = _solution(tt)[1 + _state.molefrac_offset[species.first]]
            push!(X[i], molde_fraction)
            i += 1
        end
    end

    return t, X
end

function get_energy(N)
    R = size(_solution)[1]
    t = []
    e = []

    for i in 1:R
        push!(e, [])
    end

    for tt in range(0, _tmax, N)
        push!(t, tt * t_tilde)

        push!(e[1], _solution(tt)[1] / _state.Tfrac)

        for i in 2:R
            push!(e[i], _solution(tt)[i] )
        end
    end

    return t, e
end

function get_T(N, species_name)
    Nvibmode = length(_state.species[species_name].vibmodes)
    t = []
    T = []

    for i in 1:(Nvibmode + 1)
        push!(T, [])
    end

    for tt in range(0, _tmax, N)
        push!(t, tt)
        push!(T[1], _solution(tt)[1])

        for i in 1:Nvibmode
            mole_fraciont = _solution(tt)[1 + _state.molefrac_offset[species_name]]
            theta = _state.species[species_name].vibmodes[i].theta
            degen = _state.species[species_name].vibmodes[i].degen
            f(x, p = (1, 1)) = mole_fraciont * calc_evib_kb(x, p) - _solution(tt)[i + _state.evib_offset[species_name]]
            Z = ZeroProblem(f, 1000)
            Tvib = solve(Z, Order1(), p=(theta, degen))

            push!(T[i + 1], Tvib)
        end
    end

    return t * t_tilde, T
end

function solve!(tmax)
    global _tmax = tmax / t_tilde
    problem = setup_problem!(_state, _tmax)
    global _solution =  solve(problem, alg_hints = [:stiff])
end

function initialize!(;verbose::Bool = true)
    global _state = State()
    global _verbose = verbose
    global _solution = nothing
    global _tmax = 0.0
end

function  set_verbosity!(verbose)
    global _verbose = verbose
end

function set_nrho!(nrho)
    _state.nrho = nrho
    if _verbose == true
        println("set nrho to: " * string(_state.nrho))
    end
end

function set_T!(T)
    _state.T = T
    _state.Tfrac = 1.5

    for species in _state.species
        _state.Tfrac += 0.5 * _state.mole_fractions[species.first] * species.second.dof_rot
    end

    _state.Tfrac = 1.0 / _state.Tfrac

    if _verbose == true
        println("set T to: " * string(_state.T))
    end
end

function set_Tvib!(species_name, Tvib)
    if typeof(Tvib) <: Number
        for i in eachindex(_state.Tvib[species_name])
            _state.Tvib[species_name][i] = Tvib
        end
    elseif typeof(Tvib) <: Array && size(Tvib) == size(_state.Tvib[species_name])
        for i in eachindex(_state.Tvib[species_name])
            _state.Tvib[species_name][i] = Tvib[i]
        end
    else
        error("wrong Tvib format")
    end

    if _verbose == true
        println("set Tvib of[" * species_name * "] to: " * string(_state.Tvib[species_name]))
    end 
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
