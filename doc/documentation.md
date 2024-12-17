# Energy Relaxation
The initial energy of the system is
$$
E_{tot} = E_{kin} + E_{rot} + E_{vib}
$$
Assuming that all species have the same translational and rotational temperature we can write
$$
E_{kin} = V n_{tot} \frac{3}{2} k_B T \\
E_{rot} = V k_B T \sum_s n_{s} \frac{\xi_{rot, s}}{2} \\
E_{vib} = V k_B \sum_s \left( n_{s} \sum_m g_{s, m} \frac{\theta_{s,m}}{\exp(\theta_{s,m}/T_{s,m}) - 1} \right) \\
$$
Dividng the equation by $V$ and using $\frac{n_s}{n_{tot}} = X_s$ we get
$$
e_{kin} = \frac{3}{2} k_B T \\
e_{rot} =  k_B T \sum_s X_{s} \frac{\xi_{rot, s}}{2} \\
e_{vib} = k_B \sum_s \left( X_{s} \sum_m g_{s,m} \frac{\theta_{s,m}}{\exp(\theta_{s,m}/T_{s, m}) - 1} \right) \\
$$
Since the rotational and vibrational energy have the same temperature we can write
$$
e = e_{kin} + e_{rot} = \frac{1}{2} k_B T  \left(3 +  \sum_s X_{s} \xi_s \right)
$$
with the curret temperature being calculated as
$$
T(t) = \frac{2 e(t)}{k_B T  \left(3 +  \sum_s X_{s} \xi_s \right)}.
$$
We can now calculate the change of the energy of each vibrational mode as
$$
\frac{\partial e_{m, s}}{\partial t} = \frac{e_{vib, m, s}(T) - e_{vib, m, s}(T_{m, s})}{\tau_{m, s}}
$$
where the energy per mode is
$$
e_{vib, m, s}(T) = g_{s,m} \frac{\theta_{s,m}}{\exp(\theta_{s,m}/T) - 1}.
$$
Thus the change of $e$ is
$$
\frac{\partial e}{\partial t} = -\sum_s X_s \sum_m \frac{\partial e_{m, s}}{\partial t}
$$