# Exact prescribed-flow density and adiabatic species reference

For axial v=V sin(kz), inverse characteristics satisfy
tan(k xi/2)=exp(-k V t) tan(kz/2). Their Jacobian gives
n(z,t)=n_initial(xi) dxi/dz while each specific log entropy is s_initial(xi).
The end time |k V t|=ln 2 produces a factor-two Jacobian at a stagnation point.
A separate uniform translation lasts a full axial transit. Both pressures follow
their own gamma and actual evolving density; no heat is supplied.

The flow is prescribed, with whatever mechanical forcing would maintain it left
outside this fixture. This is an exact continuity/adiabatic reference, not a full
two-fluid momentum or magnetic solution. `momentum` samples only the CURRENT
numerical density times the prescribed velocity and integrates true cell moments;
it never injects analytic evolving density or entropy into the numerical update.
The analytic functions support complex-step PDE verification at interior points.
