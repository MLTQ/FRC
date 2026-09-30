# Independent physical coupling fixtures

Uses the existing analytic cylindrical magnetic/current reference with a
radially and axially varying positive density. The regular fluid state is
initialized from its prescribed particle flux and an independent smooth
rotation; both log entropies are constant. Only the actual fluid continuity
rate supplies changing-density inertia. Analytic n_t remains an error target.

`force` is an independent smooth physical force density with all components.
Its near-axis transverse components are proportional to r. Its exact regular
momentum rate is (f_r/r,f_theta/r,f_z)/M, with inactive axis slots omitted.
This is a forced component verification, not a plasma formation trajectory.
