using DifferentialEquations
using Roots

function calc_coll_number(T::Float64, spec_data::VHS)
    spec_data.C1 * exp(spec_data.C2 * T^(-1.0/3.0)) * T^(-spec_data.omega)
end

function f(u, state, t)
    #TODO: calculate molefractions needed for calculation of Tfrac and evib

    T = u[1]
    du = zeros(length(u))
    Tfrac = calc_Tfrac(u, state)
    nrho = 0.0

    for species in state.species
        nrho += u[1 + state.nrho_offset[species.first]]
    end
    
    for reaction in state.reactions
        k = reaction.A * T^reaction.B * exp(-reaction.Ea/(kb*T))
        nu = k

        for reactant in reaction.reactants
            nrho_reactant = u[1 + state.nrho_offset[reactant]]
            nu *= nrho_reactant
        end

        nu *= t_tilde
        du[1] += nu * Tfrac * reaction.DeltaE / kb 

        for species_name in keys(reaction.stochio_coeff)
            du[1 + state.nrho_offset[species_name]] += reaction.stochio_coeff[species_name] * nu
        end
    end

    for species in state.species
        nrho_spec = u[1 + state.nrho_offset[species.first]] 
        mole_frac = nrho_spec / nrho
        nu = calc_coll_freq(species.second, nrho, T) * t_tilde
        eeq = calc_evib(T, species.second, mole_frac) / kb
        N_vibmodes = length(species.second.vibmodes)
        evib_offset = state.evib_offset[species.first]

        for v in 1:N_vibmodes
            vibmode = species.second.vibmodes[v]
            Z = state.constant_relax_mode == true ? vibmode.Z : calc_coll_number(u[1], species.second.vhs)
            tau = Z / nu
            de = (eeq[v] - u[v + evib_offset]) / tau
            du[1] -= de * Tfrac
            du[v + evib_offset] = de + du[1 + state.nrho_offset[species.first]] * u[v + evib_offset] / state.nrho
        end
    end

    return du
end

function calc_Tfrac(state)
    Tfrac = 1.5

    for species in state.species
        Tfrac += 0.5 * state.mole_fractions[species.first] * species.second.dof_rot
    end

    return 1.0 / Tfrac
end

function calc_Tfrac(u, state)
    Tfrac = 1.5
    nrho = 0.0

    for species in state.species
        nrho += u[1 + state.nrho_offset[species.first]]
    end

    for species in state.species
        mole_fraction = u[1 + state.nrho_offset[species.first]] / nrho
        Tfrac += 0.5 * mole_fraction * species.second.dof_rot
    end

    return 1.0 / Tfrac
end

function setup_problem!(state, tmax)
    u0::Vector{Float64} = []
    tspan = (0, tmax)
    evib_offset::Integer = 1
    nrho_offset::Integer = 1

    push!(u0, state.T)

    for species in state.species
        mole_fraction = state.mole_fractions[species.first]
        evib = calc_evib(state.Tvib[species.first], species.second, mole_fraction) / kb
        state.evib_offset[species.first] = evib_offset
        evib_offset += length(evib)
    
        for e in evib
            push!(u0, e)
        end
    end

    nrho_offset = evib_offset

    for species in state.species
        state.nrho_offset[species.first] = nrho_offset
        nrho_offset += 1
        push!(u0, state.mole_fractions[species.first] * state.nrho)
    end

    return ODEProblem(f, u0, tspan, state);
end

function calc_coll_freq(species, nrho, temp)
    return 4.0 * species.vhs.dref^2 * nrho * sqrt(pi * kb * species.vhs.Tref / species.mass) * (temp/species.vhs.Tref)^(1.0 - species.vhs.omega)
end

function calc_evib_kb(T, p)
    p[2] * p[1] / (exp(p[1]/T) - 1.0)
end

function calc_evib(Tvib, species, mole_frac)
    evib = []

    for i in eachindex(species.vibmodes)
        if typeof(Tvib) <: Number
            frac = species.vibmodes[i].theta / Tvib
            e = max(0.0, species.vibmodes[i].degen * kb * frac / (exp(frac) - 1.0) * Tvib)
            push!(evib, e)
        else
            frac = species.vibmodes[i].theta / Tvib[i]
            e = max(0.0, species.vibmodes[i].degen * kb * frac / (exp(frac) - 1.0) * Tvib[i])
            push!(evib, e)
        end
    end

    return mole_frac * evib
end


function calc_etot(state)
    Tfrac = calc_Tfrac(state)
    ekin = state.T / Tfrac * kb
    evib = []

    for species in state.species
        Tvib = state.Tvib[species.first]
        mole_fraction = state.mole_fractions[species.first]
        evib = vcat(evib, calc_evib(Tvib, species.second, mole_fraction))
    end

    return ekin + sum(evib)
end

function calc_etot(T, state)
    Tfrac = calc_Tfrac(state)
    ekin = T / Tfrac * kb
    evib = []

    for species in state.species
        Tvib = T
        mole_fraction = state.mole_fractions[species.first]
        evib = vcat(evib, calc_evib(Tvib, species.second, mole_fraction))
    end

    return ekin + sum(evib)
end