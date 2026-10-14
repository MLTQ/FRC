# Evolved wall current with physical carrier inertia

`WallCurrentField` carries the usual harmonic magnetic state plus y(z), the
azimuthal wall current. `WallCurrentModel` samples J=J0(B)+chi*y everywhere,
including kinetic energy, density forcing, Lorentz force and transport.
J0 uses the unchanged interior weak Ampere operator. The extra mode is L2
orthogonal to that operator's tests and is determined by the extra Ohm equation.
It carries finite current kinetic energy; it is not set by eta alone.

For E_ind restricted to the original zero-wall space, the full weak equations are

    M0 E0 = I00 J0_t + C y_t + f0,
    0 = C.T J0_t + D y_t + fchi,
    B_t=F0 E0; J0_t=A0 B_t.

Here C=I0chi, D=Ichichi, I is the actual mu/(n e²) physical current inertia;
f includes all physical forcing minus the exact external inductive lift.
Eliminate y_t=-D^-1(fchi+C.T J0_t) to use the original induction solver with
Ieff=I00-C D^-1 C.T and feff=f0-C D^-1 fchi. Reconstruction restores the full
current rate; the plasma assembly checks the complete unreduced Ohm equation.
D is positive because physical density and masses are positive. No inertia
is dropped. Only the geometry part of cached operators may be reused.

Orthogonality is exact mathematically, measured to quadrature/roundoff precision
in the unreduced equations. For uniform n/eta, the wall component relaxes on the
physical electron-current time, with the actual projected inductive forcing.
The old magnetic solution is unchanged to linear order while the missing current
component contributes consistently to momentum, heating and source work.
No exact continuum wall-Poynting equality is claimed at finite resolution.

Four complete24µs pulses, including spatial/time refinement and evolving
resistance, are recorded in docs/WALL_CURRENT_RESULTS.md. The bounded driven
component is complete with stated startup/current/model limitations. Vacuum,
free interfaces, physical coils/ends and formation remain unimplemented here.
