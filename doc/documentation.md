# Energy Relaxation
The initial energy of the system is
$$
E_{tot} = E_{kin} + E_{rot} + E_{vib}
$$
Assuming that all species have the same translational and rotational temperature we can write
$$
E_{kin} = V \frac{3}{2} k_B T \sum_s n_s \\
E_{rot} = V k_B T \sum_s n_{s} \frac{\xi_{rot, s}}{2} \\
E_{vib} = V k_B \sum_s \left( n_{s} \sum_m g_{s, m} \frac{\theta_{s,m}}{\exp(\theta_{s,m}/T_{s,m}) - 1} \right) \\
$$
Dividng the equation by $V$ and using $\sum_s n_s = n_t$ we get
$$
e_{kin} = \frac{3}{2} k_B T n_t \\
e_{rot} =  k_B T \sum_s n_{s} \frac{\xi_{rot, s}}{2} \\
e_{vib} = k_B \sum_s \left( n_{s} \sum_m g_{s,m} \frac{\theta_{s,m}}{\exp(\theta_{s,m}/T_{s, m}) - 1} \right) \\
$$
Since the rotational and vibrational energy have the same temperature we can write
$$
e = e_{kin} + e_{rot} = \frac{1}{2} k_B T  \left(3 n_t +  \sum_s n_{s} \xi_s \right)
$$
with the curret temperature being calculated as
$$
T(t) = \frac{2 e(t)}{k_B T  \left(3 n_t +  \sum_s n_{s} \xi_s \right)}.
$$

## Vibrational  relaxation
The relaxation of the vibrational energy of a species $s$ can be written as
$$
\frac{\partial e_{vib, s}}{\partial t} = \frac{\partial n_s}{\partial t} e_{m, s} + n_s \frac{\partial e_{m, s}}{\partial t}
$$
We can now calculate the change of the energy of each vibrational mode as
$$
\frac{\partial e_{m, s}}{\partial t} = \frac{e_{vib, m, s}(T) - e_{vib, m, s}(T_{m, s})}{\tau_{m, s}}
$$
where the energy per mode is
$$
e_{vib, m, s}(T) = g_{s,m} \frac{\theta_{s,m}}{\exp(\theta_{s,m}/T) - 1}.
$$
Thus the change of the $e$ is
$$
\frac{\partial e}{\partial t} = -\sum_s \frac{\partial e_{vib, s}}{\partial t} + \Delta E_R
$$
where $\Delta E_R$ is the energy from chemical reactions.