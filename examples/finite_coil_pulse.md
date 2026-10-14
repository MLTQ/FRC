# Finite-coil plasma pulse

Runs the complete moving finite-tube/current/species equations for a full24us,
with a12us prescribed sin-squared circular-coil current and12us relaxation.
R=.075m,L=1m; fully ionized Ar+,n=1e20/m3,Te/Ti=5/1eV. Two radius.12m filaments
at z=+/-.15m give100uT on axis at center. These are model inputs, not a hardware
or gas-fill recommendation. Plasma has impermeable fixed side/end boundaries.
Magnetic vacuum extends around all three interfaces to infinity; no exterior
fluid density/conductivity or artificial atmosphere exists.

The first physical comparison sets eta=1 Ohm m to expose the independent
weak-conductivity limit Jtheta=-Idot*Aunit/eta. Aunit is calculated by direct
azimuthal line integration, independently of production elliptic loop fields.
The expected integrated heat is integral Aunit²/eta dV * pi² Ipeak²/(2 Tp).
Actual finite induction, physical inertia, motion and species exchange remain.
The asymptotic comparison has physical errors as well as spatial/time errors;
eta=1 is not a plasma transport closure. Spitzer is an optional physical control.

Every accepted state/stage is retained, with frozen source and run summary.
Native thermal changes use log1p/expm1 EOS differences to resolve the pulse
above the large thermal background; no state/energy is repaired. Energy includes
full induced magnetic field and source/plasma mutual term; only the prescribed
filament self-energy baseline is excluded. Work comes from accepted Gauss stages.
Reciprocal current linkage is compared to reconstructed vacuum flux at coils.
Complete trajectory comparisons and convergence are the acceptance evidence.
No unit suite or additional verification framework. Own4GiB/256MiB guards.

Vacuum linkage samples paired (r,z) coil positions, using the returned flux
component. The initial attempt stopped in this observation setup before stepping;
its frozen source remains in finite_coil_r8_z16_t80_eta1.

For Spitzer runs the fixed-eta weak-conductivity reference is explicitly marked
inapplicable; current-error and reference-heat-normalized energy fields are null.
Actual resistance, collision/drift review diagnostics and source work are saved.

--peak-field sets positive coil-only center peak in tesla (default100uT). The
10uT Spitzer controls use the same equations, duration and residual bounds, with
energy/current scales following that actual source amplitude. The100uT case
crosses weak-field/drift review thresholds and remains explicitly conditional.
