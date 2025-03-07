using DifferentialEquations
using Roots

function calc_coll_number(T::Float64, spec_data::VHS)
    spec_data.C1 * exp(spec_data.C2 * T^(-1.0/3.0)) * T^(-spec_data.omega)
end

function f(u, state, t)
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
        eeq = calc_evib(T, species.second, nrho_spec) / kb
        N_vibmodes = length(species.second.vibmodes)
        evib_offset = state.evib_offset[species.first]
        tau_invers = 0.0

        if state.constant_relax_mode == false
            for spec_coll in state.species
                Zvib1 = calc_coll_number(T, species.second.vhs)
                Zvib2 = calc_coll_number(T, spec_coll.second.vhs)
                nu = calc_coll_freq_mix(species.second, spec_coll.second, u[1 + state.nrho_offset[spec_coll.first]], T) * t_tilde
                
                tau_invers += 1.0 / (0.5*(Zvib1 + Zvib2) / nu)
            end
        end

        for v in 1:N_vibmodes
            if state.constant_relax_mode == true
                nu = calc_coll_freq(species.second, nrho, T) * t_tilde
                tau_invers = nu / species.second.vibmodes[v].Z
            end

            de = (eeq[v] - u[v + evib_offset]) * tau_invers + du[1 + state.nrho_offset[species.first]] * u[v + evib_offset] / nrho_spec
            du[1] -= de * Tfrac
            du[v + evib_offset] = de 
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
    rot_frac = 0.0
    nrho = 0.0

    for species in state.species
        nrho += u[1 + state.nrho_offset[species.first]]
        rot_frac += species.second.dof_rot * u[1 + state.nrho_offset[species.first]]
    end

    return 2.0 / (3*nrho + rot_frac)
end

function setup_problem!(state, tmax)
    u0::Vector{Float64} = []
    tspan = (0, tmax)
    evib_offset::Integer = 1
    nrho_offset::Integer = 1

    push!(u0, state.T)

    for species in state.species
        nrho_spec = state.mole_fractions[species.first] * state.nrho
        evib = calc_evib(state.Tvib[species.first], species.second, nrho_spec) / kb
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

function calc_coll_freq_mix(species1, species2, nrho2, T)
    dref = 0.5 * (species1.vhs.dref + species2.vhs.dref)
    tref = 0.5 * (species1.vhs.Tref + species2.vhs.Tref)
    omega = 0.5 * (species1.vhs.omega + species2.vhs.omega)
    mr = species1.mass * species2.mass / (species1.mass + species2.mass)

    return 2.0 * sqrt(pi) * dref^2 * nrho2 * sqrt(2.0 * kb * tref / mr) * (T/tref)^(1.0 - omega)

    #const double mu0 = 2. * sqrt(MY_PI) * pow(d_ref, 2.) * nrho_species_[jspecies] * sqrt(2. * update->boltz * tref / mr) * pow(temp_ / tref, 1.0 - omega);
end

function calc_evib_kb(T, p)
    p[2] * p[1] / (exp(p[1]/T) - 1.0)
end

function calc_evib(Tvib, species, nrho_spec)
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

    return nrho_spec * evib
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