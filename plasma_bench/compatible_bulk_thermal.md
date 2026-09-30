# Bulk entropy, actual continuity and pressure work

For the physical log state s and the actual common density n=R n_cell,
M_n s_t=-N.T W Gamma dot grad(Ns), where Gamma is the divergence-commuting
lift of the single conservative continuity flux. Thus constant s is preserved
exactly and integral(n Ns) has zero rate up to solve/quadrature/roundoff errors,
even for compression and spatially changing n. No independent nodal continuity
equation or density projection is introduced.

Thermal energy differentiates as integral epsilon(N s_t+gamma R n_t/n).
Let h=M_n^-1 N.T W epsilon and g_n=R.T W gamma epsilon/n. The exact discrete
face-work load is D_continuity.T g_n minus L_flux.T W(h grad s). Its dot product
with the supplied faces equals thermal energy change. `face_work` is a generalized
coefficient in the stated qx/qz units; it is NOT a conservative momentum force.
The eventual momentum update must pair its kinetic-energy gradient with this
work while retaining local stress and momentum conservation. That step remains
explicitly required; no pressure heating is fitted or total energy reset.

This class uses the same scalar projection as current-driven pressure work.
Those species current terms can be added separately, with beta/-alpha factors.
Bulk compression and species transport are not a complete moving integrator,
shock entropy law or plasma-vacuum treatment. Failures leave caller state intact.
