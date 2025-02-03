module ChemicalKinetics

export set_nrho!
export set_T!
export set_Tvib!
export set_Zvib!
export set_molefrac!
export add_species!
export set_relax_mode!
export print_state
export initialize!
export execute!
export get_energies
export get_T
export get_Tvib
export get_molefrac
export get_nrho
export add_reactions!
export write2csv
export write2netCDF
export read_SPARTA_log

include("Gas.jl")
include("Reader.jl")
include("Solver.jl")
include("Constants.jl")
include("Writer.jl")

_state::State = State()
_verbose::Bool = true
_solution = nothing
_tmax::Number = 0.0

function _print(str)
    if _verbose == true
        println(str)
    end
end

function add_reactions!(file_name)
    reactions = read_reactions(file_name)

    for reaction in reactions
        push!(_state.reactions, reaction)
    end

    _print(string(length(reactions)) * " reactions added")
end

function get_nrho(N)
    t = []
    nrho = []

    for i in 1:length(_state.species)
        push!(nrho, [])
    end

    for tt in range(0, _tmax, N)
        push!(t, tt * t_tilde)
        i = 1

        for species in _state.species
            push!(nrho[i], _solution(tt)[1 + _state.nrho_offset[species.first]])
            i += 1
        end
    end

    return t, nrho
end

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
            molde_fraction = _solution(tt)[1 + _state.nrho_offset[species.first]]
            push!(X[i], molde_fraction)
            i += 1
        end
    end

    return t, X
end

function get_energies(N)
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

function get_T(N)
    t = []
    T = []

    for tt in range(0, _tmax, N)
        push!(t, tt)
        push!(T, _solution(tt)[1])
    end

    return t * t_tilde, T
end

function get_Tvib(N, species_name)
    Nvibmode = length(_state.species[species_name].vibmodes)
    
    t = []
    T = []

    for i in 1:(Nvibmode)
        push!(T, [])
    end

    for tt in range(0, _tmax, N)
        push!(t, tt)

        for i in 1:Nvibmode
            nrho_spec = _solution(tt)[1 + _state.nrho_offset[species_name]]
            theta = _state.species[species_name].vibmodes[i].theta
            degen = _state.species[species_name].vibmodes[i].degen
            f(x, p = (1, 1)) = calc_evib_kb(x, p) - _solution(tt)[i + _state.evib_offset[species_name]] / nrho_spec
            Z = ZeroProblem(f, 1000)
            Tvib = solve(Z, Order1(), p=(theta, degen))

            push!(T[i], Tvib)
        end
    end

    return t * t_tilde, T
end

function execute!(tmax)
    # check reations
    reactions = []
    N_pre = length(_state.reactions)

    for reaction in _state.reactions
        good = true

        for species_name in keys(reaction.stochio_coeff)
            if !(species_name in keys(_state.species))
                good = false
                break
            end
        end

        if good == true
            push!(reactions, reaction)
        end
    end

    _state.reactions = reactions

    _print(string(N_pre - length(_state.reactions)) * " reactions removed due to missing species")

    # check species for polyatomic speices in combnation with relax mode
    for species in _state.species
        if length(species.second.vibmodes) > 1
            _print("WARNING: variable collision numbers are not implemented for polyatomic species")
            set_relax_mode!("constant")
            break
        end
    end

    global _tmax = tmax / t_tilde
    problem = setup_problem!(_state, _tmax)
    global _solution =  solve(problem, alg_hints = [:stiff])

    _print("finished executing for tmax: [" * string(tmax) * "]")
end

function initialize!(;verbose::Bool = true)
    global _state = State()
    global _verbose = verbose
    global _solution = nothing
    global _tmax = 0.0

    _print("CHEMICAL KINETICS")
    _print("=================")
end

function  set_verbosity!(verbose)
    global _verbose = verbose
end

function set_nrho!(nrho)
    _state.nrho = nrho
    _print("set nrho to: " * string(_state.nrho))
end

function set_T!(T)
    _state.T = T

        _print("set T to: " * string(_state.T))

    for species in _state.species
        set_Tvib!(species.first, T)
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

    _print("set Tvib of[" * species_name * "] to: " * string(_state.Tvib[species_name]))
end

function set_Zvib!(species_name, Zvib)
    species = _state.species[species_name]

    if typeof(Zvib) <: Number
        for mode in species.vibmodes
            mode.Z = Zvib
        end
    elseif typeof(Zvib) <: Array && size(Zvib) == size(species.vibmodes)
        for i in eachindex(species.vibmodes)
            species.vibmodes[i].Z = Zvib[i]
        end
    else
        error("wrong Zvib format")
    end

    _print("set Zvib of[" * species_name * "] to: " * string(Zvib))
end

function set_molefrac!(species_name, mole_frac)
    _state.mole_fractions[species_name] = mole_frac
    _print("set mole fraction of [" * species_name * "] to: " * string(mole_frac))
end

function set_relax_mode!(mode::String)
    if mode == "constant"
        _state.constant_relax_mode = true
    elseif mode == "variable"
        _state.constant_relax_mode = false
    else
        error("undefined relax mode: [" * mode * "]")
    end

    _print("set relax mode to: [" * mode * "]")
end

function add_species!(file_name; mole_frac = 0.0)
    species = read_species(file_name)
    _add_species!(_state, species)
    _print("added species [" * species.name * "]")
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
