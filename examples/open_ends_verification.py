"""Open (outflow) ends: uniform magnetized axial flow against its exact downstream solution and mass outflow."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.coupled_benchmarks import aligned_meshes
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.mhd_state import from_primitive
from plasma_bench.conservative_coupled import ConservativeVacuumMHD, ConservativeParameters, ConservativeState


def run(cells=16, speed=.5, end=.4, dual_energy=False):
    mesh, grid = aligned_meshes(cells)
    sources = VacuumField(grid, (), bias=np.sqrt(MU0))      # uniform 1 T-scaled field, no current
    offset = round((grid.half_length-mesh.half_length)/mesh.dz)
    total = sources.at(0.).psi[:mesh.nr+1, offset:offset+mesh.nz+1]
    state = from_primitive(mesh, total, 1., 1., vz=speed)
    model = ConservativeVacuumMHD(mesh, sources, ConservativeState(state.conserved, np.zeros(grid.shape)),
                                  ConservativeParameters(open_ends=True, vacuum_boundary="multipole", dual_energy=dual_energy))
    model.advance(end)
    q = model._primitive(model.state, model.time)
    fast = np.sqrt(5/3*1+1)                                     # c_f = sqrt(gamma p/rho + B^2/(mu0 rho)) = sqrt(8/3)
    head = -1+(speed+fast)*model.time                           # upstream rarefaction head (lab frame)
    untouched = mesh.z > head+3*mesh.dz
    d = model.diagnostics()
    # downstream face mass flux alone, integrated: rho v pi R^2 t exactly while that end is untouched
    return dict(cells=cells, time_s=model.time, rarefaction_head_z=head,
                untouched_max_density_error=float(abs(q[:, untouched, 0]-1).max()),
                untouched_max_velocity_error=float(abs(q[:, untouched, 2]-speed).max()),
                outflow_mass_kg=d["outflow_mass_kg"], exact_outflow_mass_kg=speed*np.pi*mesh.radius**2*model.time,
                mass_relative_change_with_outflow=d["mass_relative_change"],
                fluid_energy_residual_j=d["fluid_energy_balance_residual_j"])


def main(out=Path("outputs/open_ends_verification")):
    out.mkdir(parents=True, exist_ok=True)
    rows = [run(16), run(32), run(16, dual_energy=True)]
    (out/"results.json").write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(json.dumps({k: (round(v, 12) if isinstance(v, float) else v) for k, v in r.items()}))


if __name__ == "__main__" and "--conduction" not in sys.argv:
    main(*(sys.argv[1:2] and [Path(sys.argv[1])]))


def conduction_check(cells=16, end=2e-7):
    """Open-line loss in a uniform resting plasma (psi > 0 everywhere: all lines open) against the
    independently integrated ODE dU/dt = -(U-U_w) f v_te(U)/(1.5 L), early enough that flow is negligible."""
    from scipy.integrate import solve_ivp
    from scipy.constants import e as QE, m_e as ME
    mesh, grid = aligned_meshes(cells)
    sources = VacuumField(grid, (), bias=1.)
    offset = round((grid.half_length-mesh.half_length)/mesh.dz)
    total = sources.at(0.).psi[:mesh.nr+1, offset:offset+mesh.nz+1]
    mass, n0, t0, tw, f = 3.344e-27, 1e20, 100., 2., .1
    rho, p = n0*mass, 2*n0*QE*t0
    model = ConservativeVacuumMHD(mesh, sources, ConservativeState(from_primitive(mesh, total, rho, p).conserved, np.zeros(grid.shape)),
                                  ConservativeParameters(vacuum_boundary="multipole", spitzer_ion_mass=mass,
                                                         open_line_flux_limit=f, open_line_wall_temperature_ev=tw))
    model.advance(end)
    q = model._primitive(model.state, model.time)
    out = []
    for j in (0, mesh.nz//4, mesh.nz//2):
        length = mesh.half_length-abs(mesh.z[j])+.5*mesh.dz
        uw = n0*2*QE*tw/(2/3)
        rhs = lambda t, u: [-(u[0]-uw)*f*np.sqrt(QE*u[0]*(2/3)/(2*n0*QE)/ME)/(1.5*length)]
        exact = solve_ivp(rhs, (0, model.time), [p*1.5], rtol=1e-11, atol=1e-6).y[0, -1]
        out.append(dict(z=float(mesh.z[j]), connection_m=float(length), solver_thermal_j_m3=float(q[0, j, 3]*1.5),
                        ode_thermal_j_m3=float(exact), relative_difference=float(q[0, j, 3]*1.5/exact-1)))
    d = model.diagnostics()
    return dict(cells=cells, steps=model.steps, time_s=model.time, cells_checked=out,
                open_line_loss_j=d["open_line_loss_j"], fluid_energy_residual_j=d["fluid_energy_balance_residual_j"])


if __name__ == "__main__" and "--conduction" in sys.argv:
    print(json.dumps(conduction_check(), indent=1))
