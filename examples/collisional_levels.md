# Complete resolved excitation and thermal-reservoir evolution

Accepts explicit personal local Ar I/Ar II ADF04 paths; raw provider tables and
rate arrays are not bundled, downloaded or copied into the output archive.
Input hashes identify the calculations. Ar I uses the existing current NIST
15-level energies/A coefficients with matched supplied collision strengths;
Ar II uses40 supplied LS terms, not resolved fine-structure levels.

Each fixed-charge component evolves its full populations, finite electron
thermal energy and separately produced photon energy. ne=1e20/m^3 and one
electron per modeled atom are declared inputs. For neutral Ar, background ions
provide the unevolved balancing charge; this is not a neutral-only plasma.
The source's S ionization row is explicitly excluded only for this fixed-charge
component comparison. Real moving coupling must restore consistent charge change.

Six3.5us radiation cases start at5/10/20eV for the two ions, each at350/700 steps.
Produced photon energy is never labeled escaped or net physical cooling because
absorption/re-emission and spatial transfer are absent. The neutral data's
resonance opacity can be large. Maxwellian, no-drift electrons and fixed density
are component assumptions, not inferred validity for the moving pulse.

For an independent closed limit, spontaneous radiation is disabled explicitly.
At10eV initial electron temperature, each full collision population network
relaxes for1ms at10000/20000 steps. The independent statistical-mechanical
reference solves total energy =1.5*(ne/na)*T+sum(E_i*g_i*exp(-E_i/T))/Z.
Its Boltzmann populations and temperature are compared with evolved states.
This reference uses no collision-matrix eigenvector or fitted rates.

Complete accepted histories, time-halving population/temperature differences,
particle/energy residuals, matter entropy, exact input identities and scientific
PNGs are saved. No state is clipped, normalized or reconstructed from energy.
Rejected states leave their accepted prefix and reason. No unit-test suite or
new verification framework is introduced; this evaluates a concrete missing
atomic source before any moving-plasma integration.

Run61889 completed all16 cases in13.710s; source4096b6aae2627eef2c7ee1f8ba3d4c99310305bd319f01c64ce1d0b5d22c01db.
Both output PNGs inspected. The six radiation histories have maximum10/5ns
temperature difference4.31568e-4eV and population L1 difference4.28399e-5;
relative total-energy drift <=5.89306e-13. Closed equilibrium temperatures match
the independent reference within1.05e-10eV, but transient100/50ns differences
reach0.0242344eV. Roundoff population/entropy drift is retained, not normalized.
Full results and physical limits are in docs/COLLISIONAL_LEVELS.md; no moving
coupling, charge-changing completeness or optically thick cooling is claimed.
