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