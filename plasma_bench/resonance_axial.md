# Uniform axial-source resonance transport

`axial_transport` uses translation symmetry of a stationary homogeneous cylinder
to integrate initial source position z0 exactly. Each random radial/velocity/
frequency/flight history represents every z0 in [-L/2,L/2]. It retains a live
source interval [lo,hi]; an absorption flight updates that interval by intersecting
it with [-L/2-z_relative,L/2-z_relative]. Linear flight paths attain axial extrema
at endpoints. Thus no source that crossed an end can return as a live quantum.

End losses are saved as actual z0 intervals with affine crossing times. Radial
escape removes all surviving source positions at one physical time. Every
history must end with all its original source weight accounted for. Neither
small live intervals nor long histories are discarded. Geometric interval
intersection is not a density floor, frequency cutoff or energy repair.

The detector is parallel +x through the actual aperture intersected with the
axial slit. At each emission, intersect the live z0 interval with the viewed
interval and score its length/L times exp(-kappa*distance)/(4*pi). This exact
source integration removes the rare 50 um slit-hit sampling noise. It preserves
finite cylinder ends and the complete 3D random flight and conditional thermal
frequency distribution. It is specific to a uniform axial source and homogeneous
gas, not an assumption valid for general nonuniform FRC plasma.

`resonant_decorrelation_rate` uses ORNL-TM-4368 printed p35:
nu_c = 0.904 A (g_upper/g_ground) N lambda^3 / (6 pi^2), citing Berman & Lamb,
Phys. Rev. 187, 221 (1969). It adds nu_c/(4 pi) to the Lorentz frequency HWHM.
The excited hold remains Exp(A). Conditional on hold time t, at least one
decorrelating collision occurs with probability 1-exp(-nu_c t), so integrating
over t gives Pc=nu_c/(A+nu_c), as in the report. A collision redraws the thermal
atom and Voigt rest-frequency prescription for both the actual ray and observer.
The emission time is not reset or shortened: this collision is phase/velocity
decorrelation, not quenching or spontaneous decay. Frequency-independent rates,
complete memory loss and scalar isotropic scattering are model assumptions.
The formal narrow-line rest-frequency wings remain unbounded and are reported.

Zero decorrelation retains coherent atom-frame redistribution. Packet energy
is constant h nu0; no frequency/recoil energy exchange or atom displacement is
evolved. These limitations and surface removal persist from the original model.
The source-weight sum is a distributed photon/excitation inventory budget, not
a coupled plasma energy proof. All stored exit intervals have columns:
history_id, initial_z0_low, initial_z0_high, earliest_exit_s, latest_exit_s,
side_flag (0=end,1=side). Exit time is uniform over an end interval and a point
mass for a radial interval. The first four histories retain every proposed
flight with its remaining source interval; proposed endpoints may lie beyond
the ends, while the recorded exit intervals retain the actual crossings.

| Dependent | Contract |
|---|---|
| `examples.resonance_axial_experiment` | every initial source interval exits; exact view weights, batch statistics |
| `resonance_transport.redistribute` | coherent actual/+x detunings from the same conditional absorbing atom |

Complete independent-vacuum and experimental comparisons are recorded in
`docs/RESONANCE_EXPERIMENT.md`. Source-weight residuals are <=1.33e-15. The
42mTorr measured decay is not reproduced within stated bounds; physical-model,
detector/window and narrow-line assumptions remain to be resolved explicitly.
