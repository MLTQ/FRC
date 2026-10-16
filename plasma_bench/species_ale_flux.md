# Species flux through a moving face

`moving_flux` uses the existing physical HLLC/HLL Euler Riemann flux in the local
face rest frame, then transforms its number, momentum and energy fluxes back to
the one fixed velocity frame used by the atomic state. Input face speeds are in
units of v0=sqrt(eV/m_neutral). No density, momentum or energy is changed to make
a face state admissible.

For face speed a, first subtract a from each primitive axial relative velocity.
After the Riemann solve, transform

    F_N = F_N_rest,
    F_Pz = F_Pz_rest + (m_s/m0) a F_N,
    F_E = F_E_rest + a F_Pz_rest + (m_s/m0) a² F_N/2.

Other momentum fluxes are unchanged. The result approximates F(U)-a U, including
the pressure work p v rather than incorrectly advecting pressure with the mesh.
Multiply by v0 for physical face-speed units. Reconstruction and actual cell
positivity remain the caller's responsibilities; numerical HLLC/HLL diffusion
is the same as species_euler_flux.md. The Galilean face-rest-frame construction
is the standard moving finite-volume approach; see
[Springel (2010)](https://arxiv.org/abs/0901.4107).

The moving electrostatic component calls this at interior material faces only.
Its exterior fronts use the explicitly stated vacuum Riemann condition instead
of filling exterior cells with matter.
