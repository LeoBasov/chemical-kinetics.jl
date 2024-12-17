using DifferentialEquations
using Roots

function f(u, state, t)
    T = u[1]
    du = zeros(length(u))

    for species in state.species
        mole_fraction = state.mole_fractions[species.first]
        nu = calc_coll_freq(species.second, state.nrho, T) * t_tilde
        eeq = calc_evib(T, species.second, mole_fraction) / kb
        N_vibmodes = length(species.second.vibmodes)
        evib_offset = state.evib_offset[species.first]

        for v in 1:N_vibmodes
            vibmode = species.second.vibmodes[v]
            tau = vibmode.Z / nu
            de = (eeq[v] - u[v + evib_offset]) / tau
            du[1] -= de * _state.Tfrac
            du[v + evib_offset] = de
        end
    end

    return du
end

function setup_problem!(state, tmax)
    u0::Vector{Float64} = []
    tspan = (0, tmax)
    evib_offset::Integer = 1

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

    return ODEProblem(f, u0, tspan, state);
end

function calc_coll_freq(species, nrho, temp)
    return 4.0 * species.vhs.dref^2 * nrho * sqrt(pi * kb * species.vhs.Tref / species.mass) * (temp/species.vhs.Tref)^(1.0 - species.vhs.omega)
end

function calc_evib_kb(T, p)
    p[2] * p[1] / (exp(p[1]/T) - 1.0)
end

function calc_evib(Tvib, species, mole_fraction)
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

    return mole_fraction * evib
end

function calc_etot(state)
    ekin = state.T / state.Tfrac * kb
    evib = []

    for species in state.species
        Tvib = state.Tvib[species.first]
        mole_fraction = state.mole_fractions[species.first]
        evib = vcat(evib, calc_evib(Tvib, species.second, mole_fraction))
    end

    return ekin + sum(evib)
end

function calc_etot(T, state)
    ekin = T / state.Tfrac * kb
    evib = []

    for species in state.species
        Tvib = T
        mole_fraction = state.mole_fractions[species.first]
        evib = vcat(evib, calc_evib(Tvib, species.second, mole_fraction))
    end

    return ekin + sum(evib)
end