# Complete field/current source histories

Runs the actual multispecies Lorentz/Ampere source, then the same source coupled
to the existing full reacting/colliding argon state. All accepted particle and
electric-field states, fixed frames, constant magnetic field, budgets and exact
reference histories are saved. Private atomic input paths/hashes are recorded;
raw restricted tables and rate-containing matrices are not bundled.

Two homogeneous inert limits use physical electron and singly charged argon
masses, n=1e16/m³, Te/Ti=2/0.2eV, Bz=20mT, zero initial lab velocities and
E0=100V/m. Uniform pressure has no gradient, so the momentum equations have
the exact cold-fluid oscillation solutions despite finite thermal energy.
The integration frame [200,-100,300]m/s exercises transformed work, not a lab
boost. These are verification inputs, not lab recommendations.

For E parallel to B, omega_p²=sum(n q_s²/(epsilon0 m_s)), E=E0 cos(omega_p t),
v_s=(q_s/m_s)E0 sin(omega_p t)/omega_p. For perpendicular components use
v=vx+i vy, E=Ex+i Ey. The mode condition is

    omega = sum_s omega_ps²/(omega-Omega_s), Omega_s=q_s B/m_s.

For equal singly charged ion/electron numbers it has omega=0 and
omega±=[Omega_i+Omega_e ± sqrt((Omega_i-Omega_e)²+4 omega_p²)]/2.
The associated v_s=i(q_s/m_s)E/(omega-Omega_s). Coefficients match the stated
initial E/velocities. These closed-form modes are independent of the midpoint
solver; no second numerical time integrator supplies the reference. Each
100/200/400-step history covers 5.3 periods of its fastest relevant frequency.

The two full argon histories use the existing 90% neutral/10% seed population,
nucleus density8.961429541213399e19/m³, initialTe/Ti=10/0.1eV,20eV per-nucleus
sin² heat pulse in the first microsecond, and3.5us duration. They retain the
prior neutral/ion/electron lab velocities, physical masses and lnLambda=10,
now with self-consistent E and constant Bz=20mT. This is a homogeneous source
experiment with no spatially varying field or plasma-vacuum boundary.

Requested700/1400steps are deliberate full-source timestep comparisons, not a
claim that ns steps resolve the actual fast plasma phase. Maximum omega_p*dt
and cyclotron*dt, full field/current differences and conservation errors are
recorded. Small energy drift cannot certify underresolved field dynamics.
Every failed run retains its full accepted prefix and attempted-step diagnosis;
it is not presented as a completed physical solution. No floors or corrections.

`--method exponential` changes only the local electromagnetic source integration
to the sparse exponential action; the same analytic and complete reacting cases
are retained. The original midpoint run in outputs/atomic_electromagnetic is
preserved. Its velocity-error summary mistakenly included the undefined velocity
of an empty neutral slot after a fixed-frame transform. Occupied ion/electron
errors are corrected in occupied_velocity_readout.json from the original saved
arrays; the current runner compares only those occupied species. Field/current,
temperature and budget measurements are unaffected. No scientific state changed.
