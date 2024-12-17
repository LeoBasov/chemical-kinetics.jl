mutable struct VHS
    dref::Number
    Tref::Number
    omega::Number
    alpha::Number
    Zrotinf::Number
    Tstar::Number
    C1::Number
    C2::Number

    function VHS()
        return new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    end
end

mutable struct Vibmode
    theta::Number
    degen::Int
    Z::Int

    function Vibmode()
        new(0.0, 0, 0)
    end
end

mutable struct Species
    name::String
    mass::Number
    dof_rot::Number
    dof_vib::Number
    Zrot::Int
    vhs::VHS
    vibmodes::Vector{Vibmode}

    function Species()
        new("", 1.0, 0.0, 0.0, 0, VHS(), [])
    end
end

mutable struct Reaction
    reactants::Dict{String, Number}
    products::Dict{String, Number}
    A::Number
    B::Number
    Ea::Number
    DeltaE::Number

    function Reaction()
        new(Dict(), Dict(), 0.0, 0.0, 0.0, 0.0)
    end
end

mutable struct State
    T::Number
    Tfrac::Number
    nrho::Number
    mole_fractions::Dict{String, Number}
    Tvib::Dict{String, Vector{Number}}
    species::Dict{String, Species}
    evib_offset::Dict{String, Integer}

    function State()
        return new(1.0, 1.0, 1.0, Dict(), Dict(), Dict(), Dict())
    end
end

function _add_species!(state, species)
    state.mole_fractions[species.name] = 0.0
    state.Tvib[species.name] = zeros(length(species.vibmodes))
    state.species[species.name] = species
end