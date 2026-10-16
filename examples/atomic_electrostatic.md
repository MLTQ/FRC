# Complete spatial plasma-wave histories

Exercises actual nonlinear species number/momentum/energy transport coupled to
the longitudinal Gauss field, using physical electron and argon-mass ion fluids.
The declared inert single-ion layout disables atomic reactions and collisions;
its bound-state bookkeeping is constant in the integrated energy. No chemistry
or physical chamber/preparation recommendation follows from these inputs.

The reflecting interval has L=0.1m, Nref=1e14/m³, Ti/Te=0.2/2eV, k=2pi/L and
left normal E=0. The independent linear adiabatic two-fluid density amplitudes
satisfy d²A/dt²=-K A, with

    K = [[k² ci²+omega_pi², -omega_pi²],
         [-omega_pe², k² ce²+omega_pe²]],
    cs²=gamma e Ts/ms, omega_ps²=Nref e²/(epsilon0 ms), gamma=5/3.

Choose the larger eigenvalue by its explicit quadratic formula. Its ion/electron
amplitude ratio is omega_pi²/(K_ii-omega²). Initially velocities vanish,
densities are Nref[1+a_s cos(kx)] and pressures p_s0[1+gamma a_s cos(kx)].
Exact linear densities vary as cos(omega t); v_s=a_s omega/k sin(kx)sin(omega t).
The independent field is e Nref(a_i-a_e)/(epsilon0 k) sin(kx)cos(omega t).
Cell averages include sinc(k dx/2). Nonlinear evolution differs at O(a²), so
amplitudes1e-3 and1e-4 distinguish finite-amplitude effects from numerical error.
The linear solution is not claimed exact for the finite-amplitude Euler system.

Every16/32/64-cell trajectory covers2.3 plasma-wave periods. Explicit electron
acoustic CFL and maximum omega_p*dt=0.04*16/cells both constrain each real step;
space and time refine together. No artificial mass or imposed quasi-neutrality.
The electron/ion charge difference generates E; actual face particle fluxes
give its compatible current and work. Reflecting boundaries retain pressure
stress but transfer no particles or lab energy.

All accepted states, field values, analytic arrays and measured nucleus, charge,
Gauss and particle+field energy histories are saved. Both face-lumped electric
energy and the exact integral of the cellwise-linear Gauss field are reported,
so a quadrature convention cannot conceal physical field-energy differences.
Semidiscrete Ampere/work residuals use the actual final RHS stage of each step;
they are diagnostic only. Full RK2 energy drift is measured independently.
Failures preserve every accepted prefix and attempted-step diagnosis without
repair, population clipping or a conducting exterior. This bounded longitudinal
component does not replace the formation free boundary or full magnetic solver.

The six complete 16/32/64-cell histories are preserved in
`outputs/atomic_electrostatic`; measured errors and both energy quadratures are
recorded in `docs/ATOMIC_SPATIAL_FIELDS.md`. The source archive predates the
subsequent both-empty Euler-face division bypass and retains the original code.
