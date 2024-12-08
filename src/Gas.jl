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