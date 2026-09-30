# Joint entropy-density, momentum and physical-energy flux

The explicitly new point state is (n,n vx,n vy,n vz,C_e=n s_e,C_i=n s_i).
It is neither the old cell-mean state nor the old nodal-specific-entropy state.
For smooth ideal species, C_a obeys conservative advection. Physical h_a is
internal energy per particle, p_a=(gamma_a−1)n h_a, with the original log-entropy
pressure law and physical total mass M=mi+me. H=M n|v|²/2+n(h_e+h_i).

Let bars denote arithmetic means and superscript log logarithmic means.
For transport in direction d, Gamma=n^log bar(v_d), the momentum flux is
Gamma bar(v)+bar(p)/M e_d, and each entropy-density flux is

F_Ca=bar(v_d)[n^log bar(s_a)+(gamma_a−1)(bar(n)−n^log)
              +n^log(bar(h_a)/h_a^log−1)].

This formula is derived for our two independently advected entropy densities.
It is not asserted to be an existing published two-temperature solver.
The physical H gradient is w_n=−M|v|²/2+sum((gamma_a−s_a)h_a),
w_m=M v, w_Ca=h_a. The potential is psi=p v_d. Direct use of
[s_a]=[h_a]/h_a^log−(gamma_a−1)[n]/n^log proves
[w] dot F=[psi]. A symmetric conservative flux-difference operator therefore
conserves physical H in its declared nodal integration norm as well as all
six state integrals. Its associated local energy flux is bar(w) dot F−bar(psi).
The state is not reset from that energy, and there is no global force correction.

The kinetic-preserving momentum form is motivated by
[Chandrashekar (2013)](https://arxiv.org/abs/1209.4994); the independent two-entropy
derivation and tests are ours. The discrete energy identity does not itself
prove continuum quadrature accuracy, contact preservation, positivity or shock
admissibility. No artificial dissipation, heating, ionization or atmosphere is
present. Nonpositive or unrepresentable physical values reject explicitly.
Integration into cylindrical geometry and the compatible field/current basis
remains a separate required derivation; no full plasma adoption is implied.
