# Fixed reference frame for reactive moments

`FixedFrameAtomicMoments` evolves the unchanged homogeneous reaction equations
in a fixed translating frame with velocity U. The coordinate transform is

    P'_a = P_a - m_a N_a U
    E'_a = E_a - U·P_a + m_a N_a |U|²/2.

All N_a are particles per initial nucleus, P has the inherited momentum scale,
and E is in eV/nucleus. Initial velocities are supplied in the laboratory frame;
the model subtracts U before forming extensive moments. This avoids storing
unnecessary common bulk energy when the meaningful motion is relative motion.
It does not change populations, masses, rates, pressure, temperatures or time.

`fields` and the inherited reaction generator use velocities and energies in the
integration frame. `laboratory_state` explicitly reconstructs all extensive lab
moments for output; `observe` returns lab total momentum, energy and current.
Its current transform is J_lab = J_frame + e N (sum(q_a f_a)-f_e) U,
including any evolved charge residual without enforcing neutrality. Lab energy
and momentum budgets use the actually evolved total mass and frame momentum.

The same nonlinear exponential-midpoint method runs in these coordinates.
No velocity adjustment, population floor, energy clipping or renormalization
occurs during integration. U must remain fixed; this is not an accelerating
frame or a spatial transport model. The caller records U with every output.

The reaction closure is Galilean covariant within its stated nonrelativistic,
neglected-photon-momentum approximation: scalar photon energy depends on relative
pair motion, while heavy/electron lab energy transforms with their conserved
mass and momentum. This is not a relativistic photon transformation.

Finite floating-point range remains a limitation. Poorly chosen frames or large
relative velocities can still produce unrepresentable moments/thermal energy;
the inherited positivity/empty-species rejection remains unchanged. The class
does not silently choose a frame to make a failing state pass. It is not installed
in the moving pulse; physical small-drift and missing atomic channels remain.
