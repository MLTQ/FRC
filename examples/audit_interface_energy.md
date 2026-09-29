# Instantaneous interface-energy audit

Evaluates the derived response-defect rate in `docs/INTERFACE_ENERGY_AUDIT.md`
for high, low and locally selected fluxes at common accepted shock states.
Records outward fluid power, cell magnetic-energy rate and nodal total-flux
work separately, with the actual selected limiter support. This isolates a
direct stencil change from differences between evolved trajectories.

An independent symmetric perturbation in the selected conservative/CT
direction recompletes the magnetic vacuum and differences integrated response
energy at three perturbation sizes. The source is constant in this fixture,
which the audit asserts. No state, time, energy or limiter counter is changed
by a probe. Normal advancement retains the same five output times as the
main shock measurements. Each snapshot is persisted; an interrupted file
must not be mistaken for a complete five-snapshot study. Source and audit
fingerprints are recorded separately from the main study.
The historical cell/Riemann/RK2 method is selected explicitly; this audit's
cell-field derivative is not reused for a staggered-energy state.
