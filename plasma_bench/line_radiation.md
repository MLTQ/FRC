# Resolved atomic radiation and absorption

`RadiativeLevels` loads explicitly enumerated atomic levels and spontaneous
transition probabilities, in increasing energy order. It supplies the radiative
part of a future collisional-radiative source; it does not invent collision rates
or feed the current moving pulse. `from_file` retains the exact input hash.

For every line u→l, `source(f)` returns population transfers with event rate
`R_ul = A_ul f_u` and a separate vector of produced photons per line. Populations
may sum to the neutral fraction rather than one. The caller integrates photons
as state variables. With `E_i = h c (100 E_i[cm^-1])`, the atomic energy change
plus `sum(h nu_ul R_ul)` is zero. Electron energy is not involved in spontaneous
decay. Photons are **produced**, not necessarily escaped: only an explicitly
optically thin example may identify these quantities. There is no immediate
ground-state cascade, state floor, clipping, renormalization or energy repair.

`absorption(n, k_B T, offset)` uses the Einstein relations:

    kappa_nu = A_ul lambda²/(8 pi) [(g_u/g_l)n_l - n_u] phi(nu-nu0)
    integral phi dnu = 1
    sigma_nu = nu0/c sqrt(k_B T_atom/m_atom)
    gamma_nu = (sum_j A_uj + sum_j A_lj)/(4 pi).

The Voigt profile combines thermal Doppler broadening and natural lifetime
broadening, including the lower state's radiative lifetime. Frequency uses Hz,
not angular frequency. Offset may broadcast to (samples, lines). Negative net
absorption denotes population inversion and is not clipped. Exactly zero density
gives exactly zero absorption. This is a narrow-line coefficient; no pressure,
Stark, Zeeman, velocity-gradient or collisional broadening is modeled.

The Ar I record covers ground, four 4s and ten 4p levels, and all 31 NIST-reported
lines in that energy range. Missing rates are not physically demonstrated zeros.
In particular the two 4s metastables have no radiative sink in this truncated
model. Higher states, collisions, stepwise ionization, recombination branching,
quenching and radiative transfer must be added before laboratory predictions.

| Dependent | Contract |
|---|---|
| `examples.line_radiation` | level order, line photon inventory, joule energies, SI opacity |
| Future collisional-radiative source | sum populations conserved; emitted energy separate from escape |
