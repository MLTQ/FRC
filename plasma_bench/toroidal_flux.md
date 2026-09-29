# Three-component cylindrical ideal fluxes

Uses the poloidal normal field shared by CT, with vtheta=r Omega and Btheta=r g.
Angular-momentum flux is r(rho vtheta vn-Btheta Bn/mu0); g flux is
(g vn-Omega Bn). Both are conservative under annular divergence. No division by
zero occurs at the axis. Total energy, v dot B, fast wave speed and hoop stress
retain toroidal components. Axial face quadrature uses r=sqrt(<r²>) for these
quadratic regular toroidal terms; radial faces use their actual radius.

Rusanov dissipates jumps of the actual weighted conservative variables with
the full fast characteristic bound. It is numerical diffusion, not physical
viscosity/resistivity. Replacing normal B recomputes face energy at unchanged
thermal pressure. The same radius and primitive conventions are required by
`toroidal_reconstruction` and `toroidal_mhd`.
