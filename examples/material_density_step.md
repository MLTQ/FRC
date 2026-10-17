# Strong density contrast in a free neutral slab

This complete Euler comparison isolates the pressure/mass-resolution issue
suggested by the actual Ar9 post-source state in reaction_final_trace_transport:
neighboring occupied densities differ by8.39609 decades. It does not reproduce
that reacting charged trajectory. The fixture is neutral argon, physical mass
39.948u, T_left=T_right=0.2eV, n_left=1e14m^-3 and n_right=1e6m^-3. Multiplying
both densities by a common factor leaves these neutral Euler velocities and
temperatures unchanged. The external vacuum stays exactly empty.

Initial support is[-10,10]mm with a density/pressure step at zero. Run the same
precise P0-density/P1-velocity dynamics, physical mass matrix, minmod numerical
viscosity(.5,2) and half of the actual acoustic bound used by the reacting case.
All64/128/256 uniform meshes complete a prescribed duration0.08a/sqrt(e*T0/m).
The complete central rarefaction/contact/shock lies inside|x|<a/2, while the
exterior vacuum fans cannot yet reach this comparison region. All-domain energy,
momentum and every evolving support remain part of the calculation.

The independent exact solution uses the same Rankine–Hugoniot and isentropic
wave relations as the existing sod_exact fixture, with the specified stronger
initial states. Solve f_left(p*)+f_right(p*)=0 in log(p*) between the initial
pressures. For p<=p_k, f=2*c_k/(gamma-1)*[(p/p_k)^((gamma-1)/(2gamma))-1];
otherwise f=(p-p_k)*sqrt(2/[(gamma+1)*rho_k*(p+(gamma-1)/(gamma+1)*p_k)]).
The left wave is a rarefaction and the right wave a shock; velocities start at
zero. Contact velocity, density jumps and fan boundaries follow from those
relations. No numerical trajectory or fitted coefficient enters the reference.

Report central L1 density/velocity/pressure/temperature errors from16-point
quadrature split at actual cell faces and exact wave boundaries. Save all raw
states and actual timesteps, total energy and momentum, and the numerical
temperature maximum. Exact maximum temperature is recorded separately. Large
density contrast does not justify a temperature floor, density cutoff, changed
mass, clipping, shortened success or an atmosphere. This is a physical numerical
comparison for a concrete transport problem, not a generic test framework or a
formation experiment. No unit suite is invoked. PNG output shows the finest
actual profiles against the exact central solution.

The first output material_density_step preserves its complete64-cell raw history
and a terminal serialization_failure.json: a NumPy boolean prevented writing
the first case report, before128/256 began. The corrected driver stores an
ordinary boolean and repeats the complete family in a fresh output. This changes
report serialization only, not a numerical state, physical coefficient or timestep.

The corrected full family completes64/128/256 in378/1134/2926steps. Density
and pressure L1 errors decrease, while temperature error remains about0.91–0.93
in units of initial T times initial half-width. Finest central peak0.685239eV
remains below exact1.204242eV; the analytic shock/contact and low-temperature
rarefaction remain underresolved. The inspected figure and full histories are
in outputs/material_density_step_serialized. This does not reproduce the
actual370eV trace-ion spike or validate high-contrast thermal accuracy.
