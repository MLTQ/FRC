# Complete spatial electromagnetic wave runs

Run `.venv/bin/python -m examples.atomic_maxwell <fresh-output-directory>` with
the documented scientific plotting environment. It evolves actual nonlinear
species number, three momenta and energy with all Lorentz-force components and
longitudinal/transverse fields. No unit tests or imposed analytic evolution.

The interval is 0.1 m, with reflecting matter and conducting electromagnetic
ends. Physical argon-mass ions/electrons have N=1e16/m³, Ti/Te=0.2/2 eV. This
inert one-ion bookkeeping layout disables reactions/collisions explicitly; it is
a continuum verification case, not validated laboratory argon chemistry. The
integration frame [200,-100,300] m/s differs from the lab; fields/reference
velocities are always laboratory quantities. Light speed and masses are physical.

## Independent references

All cases use k=pi/L and full duration 2.3 periods, with no fit to computed data.
The vacuum standing mode has omega=c k and exactly zero matter everywhere.
The unmagnetized linear-polarization mode has

    omega² = c² k² + sum_s omega_ps²,
    Ex = E0 sin(kz) cos(omega t),
    By = -k E0/omega cos(kz) sin(omega t),
    v_sx = q_s E0/(m_s omega) sin(kz) sin(omega t).

The circular mode has constant Bz=20 mT and uses complex Ex+iEy, Bx+iBy.
Writing Omega_s=q_s Bz/m_s with its sign, the positive high-frequency branch
is the bracketed root of the CONTINUUM dispersion equation

    omega² - c² k² - sum_s omega_ps² omega/(omega-Omega_s) = 0.

Then E=E0 sin(kz) exp(-i omega t), B=k E0/omega cos(kz) exp(-i omega t), and
v_s=i q_s E0/[m_s(omega-Omega_s)] sin(kz) exp(-i omega t). These follow directly
from the linearized physical momentum/Ampere/Faraday equations, not a matrix
assembled from the numerical solver. Spatial magnetic/current averages use the
exact sinc(k dz/2); face E is pointwise. Particle cell temperatures initialize
uniformly, so unresolved intra-cell kinetic variance is an O(E0² dz²) initial
energy difference from a finite-amplitude continuum wave.

Vacuum uses E0=1000 V/m. Both plasma modes use 100 and 1000 V/m to distinguish
finite-amplitude nonlinear departures from mesh error. Nonlinear longitudinal
motion and charge separation from v cross B are retained; their generated Ez is
reported. Linear references are not exact finite-amplitude solutions. All modes
use 16/32/64 cells and physical CFL/gyro/plasma constraints; space/time refine
together. This is not a separate time-refinement study or collisionless kinetic
validity proof.

Every accepted state and field, analytic E/B/J, integrated charge/nuclei, both
field-quadrature energy budgets and source identity are saved. Failed histories
retain accepted prefixes and attempted-step diagnosis. The final figure compares
actual full wave profiles and convergence. No total-energy correction, charge
projection, density floor or numerical atmosphere. Wall momentum transfer is not
claimed zero, and physical boundary/formation adoption remains separate.

All 15 cases completed in both `outputs/atomic_maxwell` and the subsequent
`outputs/atomic_maxwell_exact_empty` run after removing a redundant both-empty
Euler-face division. All saved arrays agree. The measured dispersion, nonlinear
longitudinal response and distinct field-energy quadratures are recorded in
`docs/ATOMIC_SPATIAL_FIELDS.md`; no further replay is needed by default.
