# Narrow side-view resonance radiance

`SideView` implements the distant, parallel-ray +x viewing approximation used
in ORNL-TM-4368 equation 3.13: at each physical emission, score
`source_weight * exp(-kappa_nu * distance)/(4*pi)` per steradian. Distance is
the exact +x intersection with the cylindrical side. The circular aperture
projects to `y²+z² <= aperture_radius²`; the axial slit additionally limits
`|z| <= slit_width/2`. No packet is redirected, split, killed or removed.

Frequency in the +x view must use the SAME absorbing atom's velocity and
rest-frame detuning as the actual emitted ray. Reusing the actual outgoing
ray frequency would lose the direction/frequency correlation. The transport
supplies this conditional observer frequency. Initial emission uses a thermal
atom and Lorentz rest frequency; complete redistribution redraws both.

`record` accumulates weighted expected detector contributions in arrival-time
bins, including flight/c. Bins are diagnostics, not timesteps. Every emission
also contributes to a per-packet integral; contributions outside the plotted
time range are recorded separately. Source weights may be likelihood ratios
for an explicitly importance-sampled initial distribution; they are never
renormalized to impose a photon or energy budget. Independent batch statistics
must include the whole correlated packet history, not treat events as samples.

This is radiance per steradian before optical efficiencies, with a constant
small solid-angle approximation, no frequency-dependent window/grating response,
no polarization and no reflected wall light. It does not claim to resolve a
finite detector's full ray optics. Absolute detector count amplitude is unknown.
Late radiative events are still followed after the detector time bins end.

| Dependent | Contract |
|---|---|
| `resonance_transport.transport` | one score at each emission, conditioned +x frequency, stable packet IDs |
| `examples.resonance_experiment` | true source likelihood ratios; per-batch histogram and per-packet totals |
