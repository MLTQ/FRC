# Actual moving reaction and unequal-mesh physical limits

`run` first evolves a cold translating neutral slab through a complete40ns
radiative decay. The declared two-level mathematical atom has3eV excitation,
A=2e7/s, hydrogen mass, n=1e19/m³, half-width20mm and lab velocity(1000,-2000,500)m/s.
Ionization, recombination and electron excitation are disabled only in this
independently soluble neutral limit. The existing full atomic moment source
transfers actual particles, momentum and energy from upper to lower states.
There are no electrons, no electric field or pressure force, and exact vacuum.

Independent continuum solution: both supports translate by V_z t, velocity
stays V, upper number=N0 exp(-At), lower=N0-upper, produced photon energy
=3eV(N0-upper), and thermal energy stays exactly zero. Every20/40/80-step
history uses source-half / actual moving material-full / source-half. New
species supports/masses are created by the conservative projection. A fixed
translating coordinate frame retains the exact cold limit without subtracting
large lab kinetic energies. Photon production accumulates separately; no escape
or spatial transport claim is made. This is real moving/source integration in
a neutral limit, not proof for hot charged reacting plasmas.

Then the original physical argon electron-ion free expansion runs full20ns,
n_i=n_e=1e14/m³, half-width10mm, Ti/Te=.2/2eV, physical masses. Independent
ion/electron material meshes have32/64,64/128,128/256 cells; cosine labels and
explicit minmod stress C1=.5,C2=2 are retained. The shared Gauss field and
compatible force/work couple the independently moving blocks. Reactions are
disabled for this unchanged transport limit. Exact union integration compares
actual J/E at equal physical times; no grid is declared an exact reference.

Raw ragged states retain every block's faces, relative velocities, thermal energy,
particle amounts, species IDs and fixed frame for every accepted step. Separate
files retain actual times and budgets. A failed source projection saves its last
actual overlap state; no floor, clipping, retry or shortened requested run.

Usage: `python -m examples.reactive_material FRESH_OUTPUT` in the usual headless
environment. Source is frozen; PNG/PDF and JSON summarize actual complete or
failed histories. No unit tests or general verification framework. Full hot
charged reaction/transport integration, practical mesh evolution and all original
FRC gates remain open even if these independent limits succeed.
