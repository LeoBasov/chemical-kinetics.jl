# ChemicalReactions.jl

## Internals
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