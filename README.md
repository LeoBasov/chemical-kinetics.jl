# ChemicalReactions.jl

## Simulation setup
Simulation setup work in 6 steps
1. Initialization
1. Adding of species
1. Adding of reactions
1. Setting initial properties
1. Execution
1. Data processing

### Initialization
```
using ChemicalKinetics

initialize!()
```

### Adding of species
```
add_species!("data/NO.json", mole_frac = 0.2)
add_species!("data/N2.json", mole_frac = 0.2)
```
The molefractions are set to 0 if not specified.
They can also be set with
```
set_molefrac!(species_name, mole_frac)
```
later.

### Adding of reactions
```
add_reactions!("data/exchange.json")
```
The reactions can be spread through multiple files with `add_reactions!(file_name)` called multiple times.

### Setting initial properties
Initial number density is set with
```
set_nrho!(nrho)
```
Initial temperature is set with
```
set_T!(T)
```
the vibraitonal temperatures of all species are automatically set to `T` if not specified otherwise.

The vibraitonal temperature of each species can be set with
```
set_Tvib!(species_name, Tvib)
```
here `Tvib` can be either a scalar or a vector with the length of the number of vibrational nodes of the species.
In this manner each vibrational node can be given its own tempearture.

The initial mole fractions can be adjusted with
```
set_molefrac!(species_name, mole_frac)
```

### Execution
The simulation is run by calling
```
execute!(tmax)
```
where `tmax` is the simulation time.

### Data processing
Simulation results can be retrieved by the functions
```
get_T(N)
get_Tvib(N, species_name)
get_molefrac(N)
get_nrho(N)
get_energies(N)
```
where `N` is the number of data points retreived in the interval 0 to tmax.

## Example
```
initialize!()

add_species!("data/NO.json", mole_frac = 0.2)
add_species!("data/N2.json", mole_frac = 0.2)
add_species!("data/N.json", mole_frac = 0.2)
add_species!("data/O2.json", mole_frac = 0.2)
add_species!("data/O.json", mole_frac = 0.2)

add_reactions!("data/exchange.json")

set_T!(10000)
set_nrho!(1e23)

execute!(4e-5)
```


## Internals

### Assumptions
The simulaiton works with mixtures but it is asumed that $T_{rot} = T$ and $T_s = T$ for all species $s$.

### Used Values
The simulation internaly works with scaled values
$$
\tilde{e} = \frac{e}{k_B}
$$
$$
\tilde{n}_{\rho} = \frac{n}{n_0}
$$
$$
\tilde{t} = \frac{t}{t_q}
$$
where $k_B$ is the Boltzmann constant, $n_0$ the total particle density at the start of the simulation, and $t_q$ a magic number of $t_q = 10^{-6}$ s.