# Sparse uniform-flow fluid linearization

Builds the linearized original regular n,q,l,mz,Ce,Ci flux equations about
uniform n, pressures and axial velocity V. Radial mass divergence and pressure
gradient use the same declared centered flux/impermeable-boundary mathematics,
including the dependent transverse axis trace. Axial differentiation uses the
declared periodic centered stencil. No new spatial operator replaces evolution.

With Z=-partial_z, R the radial particle-divergence operator and G=-partial_x,
dn_t=R dq+Z dmz; dq_t=V Z dq+G dp/M; dl_t=V Z dl;
dmz_t=V R dq+2V Z dmz-V² Z dn+Z dp/M;
dCs_t=ss(R dq+Z dmz)-ss V Z dn+V Z dCs.
dp=sum (gamma_s-1)hs[(gamma_s-ss)dn+dCs].

Returns the full sparse matrix in the usual component-last flattened ordering;
the active-coordinate caller excludes inactive axis rows/columns. Independent
nonlinear centered-difference controls must verify both rest and axial flow.
