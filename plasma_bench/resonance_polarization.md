# Polarized electric-dipole scattering

For an isolated J=0 -> 1 -> 0 electric-dipole transition without a magnetic
field, the scattering amplitude is proportional to e_out dot e_in. The NIST
argon table records ground J=0 and the 104.8 nm upper J=1. This module implements
that dipole limit, not a hydrogen core/wing mixture or a measured argon
depolarization coefficient. Initial excitation alignment remains unknown.

`dipole_emission` samples p(n|d)=3[1-(d dot n)^2]/(8 pi) by rejection from the
uniform sphere, then projects d into the outgoing transverse plane. The
normalization is analytic: integral[1-(d dot n)^2]dOmega=8 pi/3. No history is
killed or reweighted. A real unit vector is one pure linearly polarized state;
averaging random vectors represents an incoherent mixture. After absorption,
the incoming polarization becomes the next emitting dipole, preserving memory.

An unpolarized incident beam has p(mu)=3(1+mu^2)/8, <mu^2>=2/5 and single-event
linear polarization (1-mu^2)/(1+mu^2). A fixed x dipole instead has second
direction moments (1/5,2/5,2/5). These are independent analytic limits.

`observer_stokes` returns angular I,Q,U for the +x view, Q positive along z:
3/(8 pi)*(dy^2+dz^2, dz^2-dy^2, 2dy dz). Frequency and transmission must use
the same absorbing atom as the actual ray. No circular polarization, magnetic
precession, quantum interference between lines or recoil is modeled. Complete
orientation reset at a resonant collision is an explicit additional assumption,
not inferred from the measured frequency-decorrelation rate.

Primary implementation reference for the Rayleigh phase matrix and consistent
observer sampling: [SKIRT's resonant line transfer documentation](https://skirt.ugent.be/root/_dev_lyman_alpha.html),
sections Scattering phase function and Polarization. Its hydrogen branching
recipe is NOT used. Atomic resonance context: House (1970),
https://doi.org/10.1016/0022-4073(70)90002-6 (J=0-1 normal triplet).

`resonance_axial` opts in with scattering='dipole'; the scalar default retains
its earlier law. `examples.resonance_channel_experiment` compares full finite
cell trajectories, vacuum limits and actual experimental channel counts.
