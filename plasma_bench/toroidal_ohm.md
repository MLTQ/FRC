# Toroidal currents, resistive E_r/E_z and the Hall field (2.5D conservative path)

Used by `ConservativeVacuumMHD` with `toroidal=True` (and `hall=True`). The
cell variable is g = B_θ/r, averaged with weight r dr, so g·∫r dr dz = ∫B_θ dr dz
is the toroidal flux of the cell. Its update is Faraday's law on each r–z cell:
dg/dt = [E_z(r_{i+1}) − E_z(r_i)]/∫r dr − [Ē_r(z_{j+1}) − Ē_r(z_j)]/(r_c dz),
with E_z on radial faces and Ē_r (averaged over dr) on axial faces. The cell
toroidal magnetic energy is ⟨r²⟩g²/(2μ0), exact for g uniform over the cell.

## Boundaries
Every gas boundary (wall, ends, inner tube) is insulating with an exact vacuum
outside it. In that vacuum r B_θ = μ0 I_link/(2π), where I_link is the current
threading the vessel.
- **No threading current (default).** The coil current is azimuthal and the
  vacuum regions connect through the axis beyond the tube ends, so
  r B_θ = 0 and B_θ = 0 on every boundary face.
- **Threading current.** `boundary` (= a beyond the walls and ends) is that
  vacuum value, for example a rod current on the axis inside an annular
  chamber (frc-ibl.9.1). Only an annular chamber can be threaded.

On the axis r B_θ = 0 by regularity.

## Currents (`currents`)
With a = ⟨r²⟩g, summation by parts gives
dW_θ/dt = −Σ E_z J_z V_radial − Σ E_r J_r V_axial exactly, for
- J_z = 2π dz (a_i − a_{i−1})/(μ0 V_radial),
- J_r = −2π dr (a_j − a_{j−1})/(μ0 V_axial),

with a = 0 beyond the axis and a = `boundary` (normally 0) beyond the walls and
ends. `geometry` holds the dual volumes:
2π r_f dr dz on interior radial faces, π⟨r²⟩₀dz on the axis face (so
J_z(0) = 2g₀/μ0, the regular-axis value), half cells on the wall, inner-wall
and end faces, and 2π r A_n on nodes. The node volumes match the solver's
poloidal Ohmic quadrature, so J_θ (stiffness current) pairs with E_θ the same way.

## Resistive field (`ohmic`)
E_z = η J_z, Ē_r = η J_r, with the node η averaged to the faces. The heat
Σ η J² V is exactly the toroidal magnetic energy lost, and is ≥ 0.
`face_heat_to_cells` gives each cell half of each interior face's heat and all
of a boundary face's heat, for the dual-energy entropy source.

## Poynting flux (`face_toroidal_field`)
The B_θ that multiplies E in S_r = −E_z B_θ/μ0 and S_z = E_r B_θ/μ0 is the
face average of a, divided by r_f or r_c. It is zero on the axis and equals
`boundary`/R (normally zero) on every gas boundary. Then, for any E, a cell's
gas receives exactly half the J·E of each interior face it touches, and all of
it on its boundary faces. A current-free
region therefore exchanges no energy with the gas, which is the property the
"ct" magnetic-energy mode relies on. Without a threading current no B_θ energy
crosses a boundary. With one, the boundary faces carry the vacuum's Poynting
flux, which the solver counts separately (`toroidal_boundary_work_j`).

## Hall field (`hall_fields`)
E_H = c J×B with c = w/(n e) per cell. Each product is formed where its normal
field lives, and moved to the other location by the adjoint of the averaging
that brought the current there:
- J_θ J_z B_r on radial faces (B_r is the CT face field);
- J_θ J_r B_z on axial faces;
- J_z J_r B_θ in cells.

Hence Σ J_θ E_θ V_node + Σ J_z E_z V_radial + Σ J_r E_r V_axial = 0 to roundoff
for any state: the discrete Hall field does no work and only exchanges energy
between the poloidal (variational) and toroidal magnetic energies. `work`
evaluates that sum, which the solver reports as `hall_work_j` (measured
~1e-26–1e-17 J).

**Near-axis accuracy.** The skew pairing with the cell energy ⟨r²⟩g² makes the
node average of J_r carry a relative error of 1/(2k²) at node row k, for
J_r ∝ r: 50% at the first row off the axis, 12.5% at the second, 5.6% at the
third. This is an O(h²/r) absolute error. Field norms still converge
(`examples/toroidal_verification.md`), but the flux on the first node row is
10–20% off in the whistler test.

## Density cutoff (`hall_coefficient`)
**DECLARED.** w = 1/(1+(ρ_c/ρ)⁴) removes the Hall term from near-vacuum cells.
It is the complement of the gap-resistivity weight. Its largest effect is
c ≤ 0.57 m_i/(ρ_c e) at ρ = 1.32 ρ_c, which bounds the whistler rate. Without
it, the 1%-of-fill gap floor (d_i 10× the fill value) would set a whistler step
about 100× below the fill-density one.
