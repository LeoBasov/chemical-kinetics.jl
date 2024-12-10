using DifferentialEquations

function f(u, p, t)
    return 0.1*u
end

function setup_problem(state, tmax)
    u0::Vector{Float64} = []
    tspan = (0, tmax)

    push!(u0, calc_ekin_rot(state.T, state.mole_fractions, state.species))

    for species in state.species
        evib = calc_evib(state.Tvib[species.first], species.second)
    
        for e in evib
            push!(u0, e)
        end
    end

    return ODEProblem(f, u0, tspan, state);
end

function calc_coll_freq(species, nrho, temp)
    return 4.0 * species.vhs.dref^2 * nrho * sqrt(pi * kb * species.vhs.Tref / species.mass) * (temp/species.vhs.Tref)^(1.0 - species.vhs.omega)
end

function calc_ekin_rot(T, mole_fractions, species)
    e = 1.5 * kb * T

    for spec in species
        e += 0.5 * kb * mole_fractions[spec.first] * T * spec.second.dof_rot
    end

    return e
end

function calc_Tkin_rtot(ekin_rot, mole_fractions, species)
    frac_dof = 0.0

    for spec in species
        frac_dof += mole_fractions[spec.first] * spec.second.dof_rot
    end 

    return 2.0 * ekin_rot / (3.0*kb + kb*frac_dof)
end

function calc_evib(Tvib, species)
    evib = []

    for i in eachindex(species.vibmodes)
        frac = species.vibmodes[i].theta / Tvib[i]
        e = max(0.0, species.vibmodes[i].degen * kb * frac / (exp(frac) - 1.0) * Tvib[i])
        push!(evib, e)
    end

    return evib
end