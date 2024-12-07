mutable struct Species
    name::String # species name
    #frac::Number # molar fraction
    #tvib::Vector{Number} # vibrational temperature
end

mutable struct State
    temp::Number # here Tkin = Trot
    nrho::Number
    species::Vector{Species} # list of species and their properties

    function State()
        new(273.75, 1.0, [])
    end
end