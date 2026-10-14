# Prescribed axial flux in the filled cylinder

The fixed zero increment psi at r=R means the affine coefficient Bbar(t)
prescribes enclosed flux Phi=pi R² Bbar, equivalently its mean axial field.
It does not prescribe Bz(R,t). This first driven experiment has no solved
exterior vacuum or finite coils. No atmosphere is introduced outside plasma.

The analytic lift has psi_ext=r² Bbar/2 and
E_ext,theta=-r Bbar_dot/2, giving loop voltage -pi R² Bbar_dot at the wall.
Ampere acts on the induced field because the uniform lift is curl free.
Faraday for the zero-wall-flux increment uses E_ind=E_physical-E_ext.
Therefore the original weak Ohm forcing becomes f-P.T W E_ext. The single
combined induction solve retains full current inertia, pressure, Hall,
advection and all existing constraints. Merely changing Bbar is insufficient.

`AxialFluxDrive` wraps each actual transport pair with this source and uses
the actual stage field in all Lorentz/Hall work. `FluxDrivenResistivePlasma`
keeps prescribed eta plus Coulomb exchange; `FluxDrivenSpitzerPlasma` uses
evolving eta plus exchange. Their thermal/Joule physics is otherwise unchanged.
Explicit source time is required and reentrant calls on one instance reject.
`AxialFluxPulse` is a C1 sin² pulse, zero before/after its duration, with SI
tesla/second fields. It represents a prescribed flux input, not a circuit.

Native energy includes the harmonic energy V Bbar²/(2 mu0); the global cross
energy is zero from the increment boundary flux. The source power is

    P_source = integral J dot E_ext dV + V Bbar Bbar_dot/mu0.

The first term drives induced-field/current/fluid energy, and the second is
the imposed harmonic field's energy change. No source work is relabeled as
heat. The original paired Joule term supplies heating. Diagnostics keep each
term and the actual total-energy-rate residual separately. Finite-step source
work uses the accepted time-integrator stage powers; physical energy is also
reconstructed independently in the example.

The returned electric coefficient vector is E_ind. Physical electric samples
require adding the exact analytic E_ext; projecting E_ext into the old
zero-wall-Etheta space would erase its nonzero boundary value. The original
weak Ohm check is equivalent after the explicit forcing shift.

Current coefficients still occupy the original electric test space and have
zero toroidal trace at the wall. The continuum voltage-driven current need
not. This is an explicit boundary approximation whose physical convergence
must be measured; the Bessel reference retains the nonzero boundary current.
No acceptance bound is relaxed to hide this issue. Impermeable radial fluid,
periodic axial and zero toroidal magnetic-wall traces remain. This is not a
free plasma edge, real ends, sheath model, or FRC formation calculation.

Three full driven runs are recorded in docs/DRIVEN_FLUX_RESULTS.md. Fixed-eta
Bz improves on8→16 cells, but space-time current error remains28.06%→21.14%
and heating7.87%→4.47% low. The evolving-eta case is a physical control
comparison. The wall-current approximation remains an unresolved solver issue.
