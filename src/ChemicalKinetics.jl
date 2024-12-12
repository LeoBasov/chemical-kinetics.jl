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

include("Gas.jl")
include("Reader.jl")
include("Solver.jl")
include("Constants.jl")

_state::State = State()
_verbose::Bool = true
_solution = nothing
_tmax::Number = 0.0

function get_energy(N)
    R = size(_solution)[1]
    t = []
    e = []

    for i in 1:R
        push!(e, [])
    end

    for tt in range(0, _tmax, N)
        push!(t, tt)

        for i in 1:R
            push!(e[i], _solution(tt)[i])
        end
    end

    return t, e
end

function get_T(N)
    R = size(_solution)[1]
    t = []
    T = []

    for i in 1:R
        push!(T, [])
    end

    for tt in range(0, _tmax, N)
        push!(t, tt)

        Tkin_Trot = calc_Tkin_rtot(_solution(tt)[1] * kb, _state.mole_fractions, _state.species)

        push!(T[1], Tkin_Trot)

        for i in 2:R
            theta = _state.species["CH4"].vibmodes[1].theta
            degen = _state.species["CH4"].vibmodes[1].degen
            f(x, p = (1, 1)) = calc_evib_kb(x, p) - _solution(tt)[i]
            Z = ZeroProblem(f, 1000)
            Tvib = solve(Z, Order1(), p=(theta, degen))

            push!(T[i], Tvib)
        end
    end

    return t, T
end

function solve!(tmax)
    global _tmax = tmax
    problem = setup_problem(_state, _tmax)
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
