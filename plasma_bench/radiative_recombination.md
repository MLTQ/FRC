# Radiative recombination and emitted energy

`RadiativeRecombination` evaluates an explicit complete charge ladder using
published power-law or Badnell fits. `RecombiningPlasma` combines it with the
existing `ImpactIonization`, without changing that source or live pulse code.
The latter returns `(df, dEe, dEphoton, diagnostics)`, energies in J per nucleus,
rates in s^-1 and density in m^-3. This component is not yet used by the moving
simulation. An explicit ionization-disabled option isolates RR reference runs.

For initial charge q=1..Z, let Rq=ne*fq*alphaq(Te). Each event transfers one
nucleus q→q-1 and removes one electron. Temperature is Ee/(1.5*mean_charge).
For a temperature-independent cross section and Maxwellian electrons,

    alpha(T) ∝ T^(-3/2) integral E*sigma(E)*exp(-E/T) dE
    mean_capture_energy = T*(3/2 + d ln(alpha)/d ln(T)).

Here T is kB*Te in joules. We analytically differentiate the *rate fit*, including
the temperature dependence of B in Badnell's formula. The moment is exact for
the implied Maxwellian rate function, but accuracy of a fit's derivative is
not guaranteed by accuracy of its values. That uncertainty remains separate
from numerical conservation. A negative inferred moment raises, never clips.

The source equations, including impact ionization Iz and prescribed heat H, are

    dfz = I(z-1)-Iz + R(z+1)-Rz
    dEe = H - sum(chi_z*Iz) - sum(mean_capture_energy_q*Rq)
    dEphoton = sum((chi_(q-1)+mean_capture_energy_q)*Rq).

Thus Ee + atomic binding inventory + escaped photons changes only by H. The
binding release is *not* returned to electron heat. Recombination can increase
the remaining electrons' temperature by preferentially removing slow electrons
while their total thermal energy still falls. No extra 3/2*T sink is added.

The ground-state model assumes every captured electron cascades immediately
to the lower ion's ground state and all emitted photons escape. This is a
zero-density, optically thin component limit, not established for the lab tube.
It omits finite cascade times/metastable storage, reabsorption, stimulated
processes, dielectronic/three-body recombination, excitation, molecular channels,
finite ion recoil, neutral motion and species momentum/current transfer.
Do not install it as a complete moving-plasma closure.

Argon data combines modern q>=3 rates and explicitly identified older q=1,2
fits. Unknown validity intervals stay unknown and their event fraction is
reported separately from known-range extrapolation. See the data README.
Number/energy state floors, clipping, normalization and repair are absent.
