# Complete physical limits for separated density and velocity meshes

Runs actual MaterialSubcells motion with four density subcells per initial
mechanical element. All original physical masses and full durations are retained.
No unit tests, fitted trajectories, atmosphere, hidden retries or energy repairs.

Three cold translating neutral radiative histories use20/40/80 steps through40ns,
the same declared3eV,A=2e7/s mathematical atom and independent exponential decay,
translation and photon-energy formulas as reactive_material.md. Eight velocity
elements and32 density subcells initially become actual new lower/upper supports.
Only this soluble neutral limit disables ionization/recombination/electron rates.

The warm charged Ar-ion slab uses n=1e12/m³,a=25mm,T=1eV[1-(X/a)²], exact fine
cell average temperature, remote E_left=-e n a/epsilon0 and exact exterior vacuum.
The independent scalar energy integral in atomic_moving.md gives the full time
to expansion scale1.5 and its speed. Final actual position and velocity errors
are measured at32/64/128 velocity elements. Pressure and electric work retain
their physical signs; gamma=5/3 and minmod numerical Q=(.5,2) are disclosed.

Electron-ion expansion uses physical Ar ion/electron masses,n_i=n_e=1e14/m³,
a=10mm,Ti/Te=.2/2eV for the full20ns at32/64/128 mechanical elements. Atomic
reactions are disabled only for this unchanged transport comparison. Exact
union quadrature compares actual J/E at equal times; no mesh is an exact answer.

All accepted states and energy/momentum histories are saved. Fixed-number runs
store constant factored amounts/reference geometry once, and every evolved
fine face, coarse velocity and normalized thermal state. Reactive runs use
write_states with their changing population scales. Raw normalized keys and
scales permit reconstruction even if individual physical SI products underflow.

Use `python -m examples.material_subcells FRESH_OUTPUT`. Sources are frozen;
figures and JSON retain failures. This measures the concrete new solver sector,
not a successful full charged reacting argon run or FRC formation.
Failed plasma prefixes are labeled with their actual time and excluded from
equal-time spatial comparisons. This postprocessing guard was added after the
complete c0d6ab56 run; all its recorded cases reached their full requested times.
