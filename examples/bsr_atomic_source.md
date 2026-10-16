# BSR ground-argon ionization in the full atomic source

Runs the existing 72-population coupled atomic model for the same three full
3.5 µs protocols: initially Ar+, 10% ionized seed, and that seed with a declared
20 eV/nucleus sin² heat pulse. N=8.961429541213399e19 m^-3 and initial Te=10 eV.
This is a homogeneous source study, not physical gas preparation or formation.

Only the ground-neutral Ar ionization coefficient changes: the former Voronov
fit is replaced by the finite-table BSR Maxwellian cross-section integral.
Existing NIST energy levels, ADAS excitation/selected Ar II ionization,
recombination, missing channels and produced-photon assumptions remain the same.
The new excitation tables are deliberately not merged: several BSR thresholds
differ from the current physical levels, and their target identities/energy
conventions need an explicit review. No process is silently counted twice.

`argon_models` and `evolve` are reused directly. The new model runs full
350/700/1400-step histories (10/5/2.5 ns) plus the matched 1400-step old-rate
control for each protocol. Comparisons use actual common saved times. Original
number, charge, energy and external-heat diagnostics remain; failures preserve
accepted prefixes and are not turned into convergence evidence.

An independent closed-form Maxwellian integral for linear sigma(E)=a+bE on a
finite interval checks normalization and quadrature orders 4/8/16. This is a
mathematical quadrature reference, not argon physics. Actual table-integral
8/16 order differences are measured separately from time differences and from
the source-table interpolation and unsupplied high-energy contribution.

The upper energy is 300 eV. No values beyond it are invented; the report gives
the maximum speed-weighted Maxwellian fraction above it over each trajectory.
That fraction is not a bound on omitted reaction rate. No validated temperature
interval is supplied by this table; the replacement range mask is false, and
the existing outside-fit-range event diagnostic therefore includes its events.
This explicitly unknown range is separate from the known energy-domain coverage.

All personal provider paths must be supplied explicitly. Raw tables and rate
coefficient arrays stay outside Git; the source closure, hashes, accepted states,
comparisons and PNG/PDF are saved. Attribution and restrictions are in
`docs/LXCAT_BSR_INPUT.md`. No unit tests or new verification framework are used.
