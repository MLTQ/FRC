"""Annular chamber (coaxial inner tube) on the conservative MHD path, against independent solutions.

1. Static equilibrium: uniform plasma in a uniform field between an inner tube and the wall must stay at rest.
2. Resistive flux penetration: a uniform applied field ramps outside a long, heavy (effectively stationary)
   plasma annulus of uniform resistivity eta. Flux reaches the vacuum inside the inner tube only by diffusing
   through the annulus. Reference: the 1D problem in x = r^2/2, independent of the 2D code,
       d psi/dt = 2 (eta/mu0) x d2psi/dx2   for x_in < x < x_w,
   with the vacuum regions linear in x (uniform Bz): psi = psi(x_in) x/x_in inside the tube and
   psi = psi_w + (psi_b-psi_w)(x-x_w)/(x_b-x_w) in the outer gap up to the fixed induced-flux boundary r_b,
   where psi_b = B(t) x_b. Continuity of Bz = d psi/dx at both walls closes the problem. It is solved by
   second-order finite differences on 4000 nodes (Radau, rtol 1e-11) and compared at the midplane.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.sparse import diags
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.geometry import Grid
from plasma_bench.electromagnetics import VacuumField, Field
from plasma_bench.mhd_state import from_primitive
from plasma_bench.conservative_coupled import ConservativeVacuumMHD, ConservativeParameters, ConservativeState

R_IN, R_W, R_B, HALF, END_GAP = .04, .10, .12, .5, .04     # [m]
ETA, B0, DB, RAMP_S = 2e-3, .05, .2, 2e-6                    # [ohm m], [T], [T], [s]


class UniformRamp(VacuumField):
    """Uniform applied Bz(t) = B0 + DB sin^2(pi t/(2 RAMP_S)), constant after the ramp; no coils."""
    def __init__(self, grid, b0=B0, db=DB, ramp=RAMP_S):
        super().__init__(grid, ())
        self.b0, self.db, self.ramp = b0, db, ramp

    def field(self, t):
        s = min(max(t, 0.)/self.ramp, 1.)
        rate = self.db*np.pi/(2*self.ramp)*np.sin(np.pi*s) if 0. < t < self.ramp else 0.
        return self.b0+self.db*np.sin(.5*np.pi*s)**2, rate

    def at(self, t):
        b, rate = self.field(t)
        psi, dpsi = b*self.r**2/2, rate*self.r**2/2
        return Field(psi, np.zeros_like(psi), np.full_like(psi, b), dpsi,
                     np.divide(-dpsi, self.r, out=np.zeros_like(dpsi), where=self.r > 0))


def meshes(nr):
    dr = (R_W-R_IN)/nr
    dz = 4*dr
    fluid = MHDMesh(R_W, HALF, nr, int(round(2*HALF/dz)), R_IN)
    grid = Grid(R_B, HALF+END_GAP, int(round(R_B/dr))+1, int(round(2*(HALF+END_GAP)/dz))+1)
    return fluid, grid


def parameters(eta):
    return ConservativeParameters(resistivity=eta, resistivity_model="uniform", vacuum_boundary="fixed",
                                  energy_quadrature="staggered", boundary_poynting="variational", time_order=3,
                                  magnetic_energy_flux="ct", dual_energy=True, dual_energy_mode="reset",
                                  dual_energy_criterion="entropy", local_trace_bounds=True, max_step_halvings=6)


def model_for(nr, sources_kind, rho, p, eta):
    fluid, grid = meshes(nr)
    sources = UniformRamp(grid) if sources_kind == "ramp" else VacuumField(grid, (), bias=.1)
    i0, j0 = round(R_IN/fluid.dr), round(END_GAP/fluid.dz)
    total = sources.at(0.).psi[i0:i0+nr+1, j0:j0+fluid.nz+1]
    state = from_primitive(fluid, total, rho, p, energy_quadrature="staggered")
    return ConservativeVacuumMHD(fluid, sources, ConservativeState(state.conserved, np.zeros(grid.shape)), parameters(eta))


def static(nr=12, end=2e-6):
    """Uniform plasma (rho 1e-6 kg/m^3, p 100 Pa) in a uniform 0.1 T field, ideal: v_A = 89 km/s."""
    model = model_for(nr, "bias", 1e-6, 100., 0.)
    model.advance(end)
    q = model._primitive(model.state, model.time)
    d = model.diagnostics()
    return dict(nr=nr, steps=model.steps, time_s=model.time, max_speed_m_s=float(np.hypot(q[..., 1], q[..., 2]).max()),
                max_density_relative_change=float(abs(q[..., 0]/1e-6-1).max()),
                max_pressure_relative_change=float(abs(q[..., 3]/100.-1).max()),
                max_induced_flux_wb=float(abs(model.state.induced_psi).max()),
                response_energy_balance_relative=d["response_energy_balance_relative"])


def reference(times, nodes=4000):
    """Independent 1D solution: psi(x) on the annulus at the requested times."""
    x_in, x_w, x_b = R_IN**2/2, R_W**2/2, R_B**2/2
    x = np.linspace(x_in, x_w, nodes)
    h = x[1]-x[0]
    d = ETA/MU0
    ramp = UniformRamp(Grid(1., 1., 5, 5))
    main, lower, upper = -2*np.ones(nodes), np.ones(nodes-1), np.ones(nodes-1)
    upper[0] = 2.; lower[-1] = 2.                     # mirrored ghost nodes
    main[0] += -2*h/x_in                              # psi_-1 = psi_1 - 2h psi_0/x_in
    main[-1] += -2*h/(x_b-x_w)                        # psi_N+1 = psi_N-1 + 2h (psi_b - psi_N)/(x_b-x_w)
    a = diags(2*d*x/h**2)@diags([lower, main, upper], [-1, 0, 1])
    a = a.tocsc()
    edge = 2*d*x[-1]/h**2*2*h/(x_b-x_w)

    def rhs(t, y):
        out = a@y
        out[-1] += edge*ramp.field(t)[0]*x_b
        return out
    sol = solve_ivp(rhs, (0., max(times)), B0*x, method="Radau", jac=a, t_eval=times, rtol=1e-11, atol=1e-14)
    return x, sol.y.T


def profile(x, psi_annulus, b):
    """Reference psi(r) on [0, R_B] from the annulus solution and the two linear-in-x vacuum regions."""
    x_in, x_w, x_b = R_IN**2/2, R_W**2/2, R_B**2/2
    spline = CubicSpline(x, psi_annulus)

    def psi(r):
        xr = np.asarray(r)**2/2
        out = spline(np.clip(xr, x_in, x_w))
        out = np.where(xr < x_in, psi_annulus[0]*xr/x_in, out)
        return np.where(xr > x_w, psi_annulus[-1]+(b*x_b-psi_annulus[-1])*(xr-x_w)/(x_b-x_w), out)
    return psi


def diffusion(nr=12, times=(1e-6, 2e-6, 3e-6, 4e-6), rho=100., p=1e3):
    """Heavy plasma (rho 100 kg/m^3, speeds < 0.03 m/s): v x B < 1e-4 eta J even after the ramp, when J is small
    (at rho = 1 kg/m^3 the residual motion was a ~1% effect at 4 us)."""
    model = model_for(nr, "ramp", rho, p, ETA)
    x, ref = reference(np.array(times))
    x2, ref2 = reference(np.array(times), nodes=2000)
    grid = model.sources.grid
    j = grid.nz//2
    rows, start = [], time.time()
    for t, psi_ref, psi_ref2 in zip(times, ref, ref2):
        model.advance(t)
        b = model.sources.field(model.time)[0]
        psi = model.total_flux()[:, j]
        exact = profile(x, psi_ref, b)(grid.r)
        induced = exact-b*grid.r**2/2
        inner_field = 2*psi[round(R_IN/grid.dr)]/R_IN**2
        d = model.diagnostics()
        q = model._primitive(model.state, model.time)
        rows.append(dict(nr=nr, time_s=model.time, steps=model.steps, applied_t=b,
                         inner_field_t=float(inner_field), inner_field_reference_t=float(psi_ref[0]/(R_IN**2/2)),
                         wall_field_t=float((psi[-1]-psi[-2])/(.5*(grid.r[-1]**2-grid.r[-2]**2))),
                         max_flux_error_relative=float(abs(psi-exact).max()/abs(induced).max()),
                         reference_self_error_relative=float(abs(profile(x2, psi_ref2, b)(grid.r)-exact).max()/abs(induced).max()),
                         max_speed_m_s=float(np.hypot(q[..., 1], q[..., 2]).max()),
                         response_energy_balance_relative=d["response_energy_balance_relative"],
                         ohmic_heat_j=d["ohmic_heat_j"], source_work_j=d["source_work_j"],
                         mass_relative_change=d["mass_relative_change"], wall_s=time.time()-start))
        print(json.dumps(rows[-1]), flush=True)
    return rows


def main(out=Path("outputs/annulus_verification"), resolutions=(12, 24)):
    out.mkdir(parents=True, exist_ok=True)
    result = dict(static=static(), diffusion={str(n): diffusion(n) for n in resolutions})
    print(json.dumps(result["static"]))
    errors = [result["diffusion"][str(n)][-1]["max_flux_error_relative"] for n in resolutions]
    result["observed_order"] = [float(np.log2(a/b)) for a, b in zip(errors[:-1], errors[1:])]
    (out/"results.json").write_text(json.dumps(result, indent=1))
    print(json.dumps(dict(errors=errors, observed_order=result["observed_order"])))


if __name__ == "__main__":
    main(*(sys.argv[1:2] and [Path(sys.argv[1])]),
         **({"resolutions": tuple(int(n) for n in sys.argv[2].split(","))} if len(sys.argv) > 2 else {}))
