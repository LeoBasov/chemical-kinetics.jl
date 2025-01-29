mutable struct Reaction
    stochio_coeff::Dict{String, Float64}
    reactants::Vector{String}
    A::Float64
    B::Float64
    Ea::Float64
    DeltaE::Float64

    function Reaction()
        new(Dict(), [], 0.0, 0.0, 0.0, 0.0)
    end
end

mutable struct VHS
    dref::Float64
    Tref::Float64
    omega::Float64
    alpha::Float64
    Zrotinf::Float64
    Tstar::Float64
    C1::Float64
    C2::Float64

    function VHS()
        return new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    end
end

mutable struct Vibmode
    theta::Float64
    degen::Int
    Z::Int

    function Vibmode()
        new(0.0, 0, 0)
    end
end

mutable struct Species
    name::String
    mass::Float64
    dof_rot::Float64
    dof_vib::Float64
    Zrot::Int
    vhs::VHS
    vibmodes::Vector{Vibmode}

    function Species()
        new("", 1.0, 0.0, 0.0, 0, VHS(), [])
    end
end

mutable struct State
    T::Float64
    nrho::Float64
    mole_fractions::Dict{String, Float64}
    Tvib::Dict{String, Vector{Float64}}
    species::Dict{String, Species}
    evib_offset::Dict{String, Integer}
    molefrac_offset::Dict{String, Integer}
    reactions::Vector{Reaction}
    constant_relax_mode::Bool 

    function State()
        return new(1.0, 1.0, Dict(), Dict(), Dict(), Dict(), Dict(), [], true)
    end
end

function _add_species!(state, species)
    state.mole_fractions[species.name] = 0.0
    state.Tvib[species.name] = zeros(length(species.vibmodes))
    state.species[species.name] = species
end