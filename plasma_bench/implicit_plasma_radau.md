# Full-equation Radau IIA candidate

Uses the same validated physical coordinate/preconditioner construction as
Gauss, with two-stage Radau IIA A=[[5/12,-1/12],[3/4,1/4]],
b=[3/4,1/4], c=[1/3,1]. Both actual nonlinear stages satisfy the original
coupled equations; the second stage is the accepted endpoint. The strict
collocation engine retains all original nonlinear/induction acceptance bounds,
physical states and atomic failures. No fitted fields or thermodynamic repairs.

R(z)=(1+z/3)/(1-2z/3+z²/6), third order and L stable. Unlike Gauss, this
time discretization intentionally damps unresolved oscillations; it is not an
energy-conserving or dissipation-free evolution method. Physical fast waves
are also attenuated when unresolved. Report wave amplitude/phase and energy
versus time step and compare full contact/mesh controls. Damping cannot be
presented as physical resistivity or evidence of physical plasma stability.
This remains a bounded filled-cylinder verification candidate until those
controls pass; it does not implement vacuum or formation.
