"""2.5D (v_theta, B_theta) and Hall verification of ConservativeVacuumMHD against independent exact solutions.

Run: python -m examples.toroidal_verification OUT CASE [nr ...]   CASE in torsional, rotor, alfvenic, whistler, rod.
All cases use the formation solver settings (staggered energy, variational boundary Poynting flux, RK3, CT
magnetic-energy flux, local trace bounds) without dual energy, in a gas cylinder with insulating walls and ends
and the exact vacuum outside. Cell averages of every reference are taken with 4x4 Gauss points per cell.

1. torsional: standing torsional Alfven wave in a uniform axial field (ideal). With B_r = 0 every field line
   is independent, so for any radial profile f with f(R) = 0
       v_theta = V f(r) cos(k s) cos(w t),   B_theta = -sqrt(mu0 rho) V f(r) sin(k s) sin(w t),
   s = z + L, k = pi/(2L), w = k v_A: B_theta vanishes at both insulating ends (free ends: v_theta antinodes).
   f = x(1-x^2)^3, x = r/R. V = 1e-3 v_A keeps the second-order (compressive) coupling below 1e-6.
2. rotor: isothermal rigid rotor in a uniform axial field, v_theta = W r, p = rho c^2,
   rho = rho0 exp(W^2 r^2/(2 c^2)): an exact steady state (centrifugal force against the pressure gradient).
3. alfvenic: a nonlinear stationary Alfvenic state with no poloidal field, v = B/sqrt(mu0 rho) along a purely
   toroidal field, B_theta = B1 f(r) cos^3(pi z/(2L)), p = P0 - B_theta^2/(2 mu0), uniform rho. Hoop stress
   cancels the centrifugal force and the total pressure is uniform, so it tests the B_theta pressure, hoop and
   centrifugal terms together, at Mach ~1 and beta ~1.
4. whistler: electron-MHD (Hall) evolution of a toroidal-field packet in a uniform B0 with immobile ions
   (ion mass 1e6 m_p: v_A = 220 m/s, w_ci/w ~ 1e-6). The exact linear solution in unbounded space is a
   Hankel-Fourier superposition with w = D |k| |k_z|, D = B0/(mu0 n e) (parallel limit w = k^2 d_i v_A):
       b_theta = b0 cos(w t),  b_r = -sgn(D) (|k_z|/|k|) b0 sin(w t),  b_z = i k_perp b_r / k_z,
   evaluated by Gauss quadrature, independent of the solver. The packet B_theta = A (r/w_r) exp(-r^2/w_r^2)
   exp(-z^2/w_z^2) cos(k0 z) stays e^-20 below its peak at every boundary over one period. The measured
   whistler frequency is the scale s on D (hence on w) that best fits the solver's fields.
5. rod: a current I(t) threads an annular chamber (toroidal_boundary = mu0 I/(2 pi), ramped as sin^2 over
   1 us). B_theta diffuses resistively into a heavy, stationary annulus from both walls and the ends. At the
   midplane of the long annulus a = r B_theta obeys the 1D problem
       da/dt = (eta/mu0) r d/dr((1/r) da/dr),   a(r_in) = a(r_w) = mu0 I(t)/(2 pi),
   solved independently by second-order finite differences on 4000 intervals (Radau, rtol 1e-11). The
   separately counted boundary B_theta Poynting work must equal the toroidal magnetic energy plus Ohmic heat.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE, m_p as MP
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
from scipy.special import j0, j1
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.geometry import Grid
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.magnetic_quadrature import selected_energy
from plasma_bench.conservative_coupled import ConservativeVacuumMHD, ConservativeParameters, ConservativeState

GAMMA = 5/3
NODES, WEIGHTS = np.polynomial.legendre.leggauss(4)


def parameters(**extra):
    settings = dict(resistivity=0., vacuum_boundary="fixed", energy_quadrature="staggered",
                    boundary_poynting="variational", time_order=3, magnetic_energy_flux="ct",
                    local_trace_bounds=True, toroidal=True)
    return ConservativeParameters(**{**settings, **extra})


def cell_average(mesh, function, weight=lambda r: r):
    """Average of function(r, z) over each cell with radial weight `weight` (4x4 Gauss)."""
    total = np.zeros((mesh.nr, mesh.nz))
    norm = np.zeros((mesh.nr, 1))
    for xr, wr in zip(NODES, WEIGHTS):
        r = (mesh.r+.5*mesh.dr*xr)[:, None]
        norm += wr*weight(r)
        for xz, wz in zip(NODES, WEIGHTS):
            z = (mesh.z+.5*mesh.dz*xz)[None, :]
            total += wr*wz*weight(r)*function(r, z)
    return total/(2*norm)


def build(nr, nz, radius, half, rho, pressure, angular, g_avg, *, bias=0., **extra):
    """Model with cell density/pressure arrays and cell L = <rho r v_theta>, g = <B_theta/r> (r dr weighted)."""
    fluid = MHDMesh(radius, half, nr, nz)
    grid = Grid(radius+3*fluid.dr, half+3*fluid.dz, nr+4, nz+7)
    sources = VacuumField(grid, (), bias=bias)
    psi = sources.at(0.).psi[:nr+1, 3:3+nz+1]
    r2 = .5*(fluid.rf[:-1]**2+fluid.rf[1:]**2)[:, None]
    rho = np.broadcast_to(rho, (nr, nz)); momentum = np.broadcast_to(angular, (nr, nz))
    g = np.broadcast_to(g_avg, (nr, nz))
    energy = (pressure/(GAMMA-1)+selected_energy(fluid, psi, "staggered")+.5*momentum**2/(rho*r2)
              +r2*g**2/(2*MU0))
    conserved = np.stack((rho, np.zeros((nr, nz)), np.zeros((nr, nz)), energy, momentum, g), axis=-1)
    model = ConservativeVacuumMHD(fluid, sources, ConservativeState(conserved, np.zeros(grid.shape)),
                                  parameters(**extra))
    return model


def relative_l1(mesh, value, reference, scale):
    return float(np.sum(np.abs(value-reference)*mesh.volumes)/(np.sum(mesh.volumes)*scale))


def energy_report(model):
    d, start = model.diagnostics(), model.initial
    toroidal0 = start["toroidal_kinetic_energy_j"]+start["toroidal_magnetic_energy_j"]
    toroidal = d["toroidal_kinetic_energy_j"]+d["toroidal_magnetic_energy_j"]
    return dict(toroidal_energy_relative_change=(toroidal-toroidal0)/toroidal0 if toroidal0 else None,
                fluid_energy_residual_j=d["fluid_energy_balance_residual_j"],
                response_energy_residual_j=d["response_energy_balance_residual_j"],
                response_energy_j=d["response_energy_j"], toroidal_magnetic_energy_j=d["toroidal_magnetic_energy_j"],
                toroidal_kinetic_energy_j=d["toroidal_kinetic_energy_j"], hall_work_j=d["hall_work_j"],
                steps=d["steps"])


# 1. Torsional Alfven standing wave -----------------------------------------------------------------------
R, HALF, RHO0, B0 = .1, .1, 1e-6, .1
P0 = B0**2/(2*MU0)       # beta = 1


def profile(r):
    x = r/R
    return x*(1-x*x)**3


def torsional(nr, periods=1.):
    va = B0/np.sqrt(MU0*RHO0)
    amplitude, k = 1e-3*va, np.pi/(2*HALF)
    w = k*va
    mesh_l = lambda r, z: RHO0*r*profile(r)*np.cos(k*(z+HALF))           # rho r v_theta mode, at cos(wt) = 1
    mesh_g = lambda r, z: profile(r)*np.sin(k*(z+HALF))/r                  # B_theta/r mode
    fluid = MHDMesh(R, HALF, nr, 2*nr)
    phi_l, phi_g = cell_average(fluid, mesh_l), cell_average(fluid, mesh_g)
    model = build(nr, 2*nr, R, HALF, RHO0, P0, amplitude*phi_l, 0., bias=B0)
    vol = fluid.volumes
    series = []

    def record(m):
        u = m.state.conserved
        a_v = np.sum(vol*u[..., 4]*phi_l)/np.sum(vol*phi_l**2)
        a_b = np.sum(vol*u[..., 5]*phi_g)/np.sum(vol*phi_g**2)
        series.append((m.time, a_v/amplitude, -a_b/(np.sqrt(MU0*RHO0)*amplitude)))
    record(model)
    start = time.time()
    model.advance(periods*2*np.pi/w, on_step=record)
    t, cv, sb = np.array(series).T
    phase = np.unwrap(np.arctan2(sb, cv))
    fit = np.polyfit(t, phase, 1)
    amp = np.hypot(cv, sb)
    u = model.state.conserved
    exact_l = amplitude*phi_l*np.cos(w*model.time)
    exact_g = -np.sqrt(MU0*RHO0)*amplitude*phi_g*np.sin(w*model.time)
    return dict(case="torsional", nr=nr, nz=2*nr, time_s=model.time, wall_s=time.time()-start,
                frequency_relative_error=float(fit[0]/w-1), amplitude_change=float(amp[-1]-1),
                angular_momentum_l1=relative_l1(fluid, u[..., 4], exact_l, np.max(np.abs(amplitude*phi_l))),
                toroidal_field_l1=relative_l1(fluid, u[..., 5], exact_g, np.sqrt(MU0*RHO0)*amplitude*np.max(np.abs(phi_g))),
                **energy_report(model))


# 2. Isothermal rigid rotor ---------------------------------------------------------------------------------
def rotor(nr, crossings=3.):
    cs2 = P0/RHO0                       # isothermal sound speed^2
    spin = .5*np.sqrt(cs2)/R            # Mach 0.5 at the wall
    density = lambda r, z: RHO0*np.exp(spin**2*r**2/(2*cs2))+0*z
    fluid = MHDMesh(R, HALF, nr, 2*nr)
    rho = cell_average(fluid, density)
    pressure = cs2*rho
    momentum = cell_average(fluid, lambda r, z: density(r, z)*spin*r*r)
    model = build(nr, 2*nr, R, HALF, rho, pressure, momentum, 0., bias=B0)
    start = time.time()
    model.advance(crossings*R/np.sqrt(GAMMA*cs2))
    q = model._primitive(model.state, model.time)
    u = model.state.conserved
    return dict(case="rotor", nr=nr, nz=2*nr, time_s=model.time, wall_s=time.time()-start,
                density_l1=relative_l1(fluid, u[..., 0], rho, RHO0),
                angular_momentum_l1=relative_l1(fluid, u[..., 4], momentum, np.max(momentum)),
                pressure_l1=relative_l1(fluid, q[..., 3], pressure, P0),
                max_poloidal_speed_over_cs=float(np.hypot(q[..., 1], q[..., 2]).max()/np.sqrt(cs2)),
                **energy_report(model))


# 3. Nonlinear stationary Alfvenic toroidal state ------------------------------------------------------------
def alfvenic(nr, crossings=2.):
    b1 = B0
    peak = (1/np.sqrt(7))*(6/7)**3   # max of x(1-x^2)^3
    btheta = lambda r, z: b1*profile(r)/peak*np.cos(np.pi*z/(2*HALF))**3
    total_pressure = 2*b1**2/(2*MU0)
    fluid = MHDMesh(R, HALF, nr, 2*nr)
    g = cell_average(fluid, lambda r, z: btheta(r, z)/r)
    momentum = cell_average(fluid, lambda r, z: np.sqrt(RHO0/MU0)*btheta(r, z)*r)
    pressure = cell_average(fluid, lambda r, z: total_pressure-btheta(r, z)**2/(2*MU0))
    model = build(nr, 2*nr, R, HALF, RHO0, pressure, momentum, g, bias=0.)
    start = time.time()
    va = b1/np.sqrt(MU0*RHO0)
    model.advance(crossings*R/va)
    q = model._primitive(model.state, model.time)
    u = model.state.conserved
    return dict(case="alfvenic", nr=nr, nz=2*nr, time_s=model.time, wall_s=time.time()-start,
                density_l1=relative_l1(fluid, u[..., 0], RHO0, RHO0),
                angular_momentum_l1=relative_l1(fluid, u[..., 4], momentum, np.max(momentum)),
                toroidal_field_l1=relative_l1(fluid, u[..., 5], g, np.max(g)),
                pressure_l1=relative_l1(fluid, q[..., 3], pressure, total_pressure),
                max_poloidal_speed_over_va=float(np.hypot(q[..., 1], q[..., 2]).max()/va),
                **energy_report(model))


# 4. Whistler (electron-MHD) packet --------------------------------------------------------------------------
W_B0, W_N, W_MASS, W_AMPLITUDE = .1, 1e20, 1e6*MP, 1e-5     # [T], [m^-3], [kg], packet scale [T]
W_LAMBDA = .04
W_K0, W_WR, W_WZ = 2*np.pi/W_LAMBDA, .04, .05
W_R, W_HALF = .2, .26
W_D = W_B0/(MU0*W_N*QE)


class WhistlerReference:
    """Exact linear EMHD field of the packet (Hankel order 1/0 in r, cosine/sine in z) by Gauss quadrature."""
    def __init__(self, count=320):
        x, w = np.polynomial.legendre.leggauss(count)
        kp_max, kz_max = 14/W_WR, W_K0+14/W_WZ
        self.kp, self.wp = .5*kp_max*(x+1), .5*kp_max*w
        self.kz, self.wz = .5*kz_max*(x+1), .5*kz_max*w
        b0 = self.kp*W_WR**3/4*np.exp(-self.kp**2*W_WR**2/4)      # Hankel-1 of (r/w_r) exp(-r^2/w_r^2)
        zhat = np.sqrt(np.pi)*W_WZ/2*(np.exp(-(self.kz-W_K0)**2*W_WZ**2/4)+np.exp(-(self.kz+W_K0)**2*W_WZ**2/4))
        self.amplitude = W_AMPLITUDE*(self.wp*self.kp*b0)[:, None]*(self.wz*zhat)[None, :]/np.pi   # (kp, kz)
        self.k = np.hypot(self.kp[:, None], self.kz[None, :])

    def fields(self, r, z, t, scale=1.):
        """(B_r, B_theta, B_z) at points r (nr,), z (nz,) for Hall coefficient scaled by `scale`."""
        phase = scale*W_D*self.k*self.kz[None, :]*t
        a_t = self.amplitude*np.cos(phase)
        a_r = -np.sign(W_D)*self.amplitude*(self.kz[None, :]/self.k)*np.sin(phase)
        a_z = -self.kp[:, None]*a_r/self.kz[None, :]
        cz, sz = np.cos(np.outer(self.kz, z)), np.sin(np.outer(self.kz, z))
        bessel1, bessel0 = j1(np.outer(r, self.kp)), j0(np.outer(r, self.kp))
        return bessel1@(a_r@cz), bessel1@(a_t@cz), bessel0@(a_z@sz)

    def flux(self, r, z, t, scale=1.):
        """Poloidal flux psi = integral of r B_z dr, using integral of r J0(k r) dr = r J1(k r)/k."""
        phase = scale*W_D*self.k*self.kz[None, :]*t
        a_r = -np.sign(W_D)*self.amplitude*(self.kz[None, :]/self.k)*np.sin(phase)
        a_z = -self.kp[:, None]*a_r/self.kz[None, :]
        return (r[:, None]*j1(np.outer(r, self.kp))/self.kp[None, :])@(a_z@np.sin(np.outer(self.kz, z)))


def gauss_points(mesh):
    r = (mesh.r[:, None]+.5*mesh.dr*NODES[None, :]).ravel()
    z = (mesh.z[:, None]+.5*mesh.dz*NODES[None, :]).ravel()
    return r, z


def toroidal_average(mesh, btheta_points):
    """<B_theta/r> with r dr weight over each cell, from B_theta at the 4x4 Gauss points: integral B_theta dr dz."""
    b = btheta_points.reshape(mesh.nr, 4, mesh.nz, 4)
    return .25*np.einsum("iajb,a,b->ij", b, WEIGHTS, WEIGHTS)*mesh.dr*mesh.dz/(mesh.radial_measure[:, None]*mesh.dz)


def relative_l2(mesh, value, reference):
    return float(np.sqrt(np.sum(mesh.volumes*(value-reference)**2)/np.sum(mesh.volumes*reference**2)))


def whistler(nr, periods=1., fraction=.5):
    """Measured at a quarter period (poloidal field, through the axis B_z) and after `periods` (B_theta)."""
    reference = WhistlerReference()
    nz = int(round(2*W_HALF/(W_R/nr)))
    fluid = MHDMesh(W_R, W_HALF, nr, nz)
    r, z = gauss_points(fluid)
    g0 = toroidal_average(fluid, reference.fields(r, z, 0.)[1])
    model = build(nr, nz, W_R, W_HALF, W_N*W_MASS, 1e3, 0., g0, bias=W_B0, hall=True, hall_ion_mass=W_MASS,
                  spitzer_ion_mass=W_MASS, hall_whistler_fraction=fraction)
    w0 = W_D*W_K0**2
    start = time.time()
    model.advance(.25*2*np.pi/w0)
    # Induced poloidal flux on all gas nodes, weighted by the node radius (the flux carried by each ring).
    psi = model.state.induced_psi[:fluid.nr+1, 3:3+fluid.nz+1]
    exact = reference.flux(fluid.rf, fluid.zf, model.time)
    weight = fluid.rf[:, None]
    row1 = float(np.linalg.norm(psi[1]-exact[1])/np.linalg.norm(exact[1]))
    quarter = dict(time_s=model.time, induced_flux_relative_l2=float(np.linalg.norm(weight*(psi-exact))/np.linalg.norm(weight*exact)),
                   first_row_flux_relative_l2=row1)
    model.advance(periods*2*np.pi/w0)
    wall = time.time()-start
    # B_theta norm: g weighted by <r^2> (the toroidal field energy), so the O(h^2/r) near-axis error of
    # the skew Hall pairing (toroidal_ohm.md) is not magnified by 1/r.
    g = model.state.conserved[..., 5]
    weight = fluid.volumes*.5*(fluid.rf[:-1]**2+fluid.rf[1:]**2)[:, None]
    norm = lambda a, b: float(np.sqrt(np.sum(weight*(a-b)**2)/np.sum(weight*b**2)))
    exact = lambda scale: toroidal_average(fluid, reference.fields(r, z, model.time, scale)[1])
    fit = minimize_scalar(lambda s: np.sum(weight*(g-exact(s))**2), bounds=(.5, 1.5),
                          method="bounded", options=dict(xatol=1e-7))
    near = fluid.r < W_WR/4
    error = weight*(g-exact(1.))**2
    return dict(case="whistler", nr=nr, nz=nz, time_s=model.time, wall_s=wall, fraction=fraction,
                omega0_rad_s=w0, toroidal_field_relative_l2=norm(g, exact(1.)),
                frequency_relative_error=float(fit.x-1),
                toroidal_field_relative_l2_after_fit=norm(g, exact(fit.x)),
                near_axis_error_share=float(error[near].sum()/error.sum()),
                toroidal_over_radius_relative_l2=relative_l2(fluid, g, exact(1.)),
                quarter_period=quarter, hall_limited_steps=model.hall_limited_steps, **energy_report(model))


# 5. Current threading an annulus ---------------------------------------------------------------------------
ROD_IN, ROD_W, ROD_HALF, ROD_ETA, ROD_K0, ROD_RAMP = .04, .10, .3, 2e-3, 2e-3, 1e-6   # m, m, m, ohm m, T m, s


def rod_boundary(t):
    """r B_theta of the vacuum: mu0 I/(2 pi), I ramping to 10 kA as sin^2 over ROD_RAMP."""
    return ROD_K0*np.sin(.5*np.pi*min(max(t, 0.)/ROD_RAMP, 1.))**2


def rod_reference(times, count=4000):
    r = np.linspace(ROD_IN, ROD_W, count+1)
    h, half = r[1]-r[0], .5*(r[:-1]+r[1:])
    d = ROD_ETA/MU0

    def rhs(t, a):
        full = np.concatenate(([rod_boundary(t)], a, [rod_boundary(t)]))
        flux = np.diff(full)/(h*half)
        return d*r[1:-1]*np.diff(flux)/h
    sol = solve_ivp(rhs, (0., max(times)), np.zeros(count-1), method="Radau", t_eval=times, rtol=1e-11, atol=1e-16)
    return r, [np.concatenate(([rod_boundary(t)], sol.y[:, k], [rod_boundary(t)])) for k, t in enumerate(times)]


def rod(nr, end_us=3.):
    dr = (ROD_W-ROD_IN)/nr
    nz = int(round(2*ROD_HALF/(4*dr)))
    fluid = MHDMesh(ROD_W, ROD_HALF, nr, nz, ROD_IN)
    grid = Grid(ROD_W+3*dr, ROD_HALF+3*fluid.dz, int(round((ROD_W+3*dr)/dr))+1, nz+7)
    sources = VacuumField(grid, ())
    i0 = int(round(ROD_IN/dr))
    rho, pressure = 1., 1e5      # heavy and stiff: displacement ~1 um over the test
    energy = pressure/(GAMMA-1)*np.ones((nr, nz))
    zero = np.zeros((nr, nz))
    conserved = np.stack((rho+zero, zero, zero, energy, zero, zero), axis=-1)
    model = ConservativeVacuumMHD(fluid, sources, ConservativeState(conserved, np.zeros(grid.shape)),
                                  parameters(resistivity=ROD_ETA, toroidal_boundary=rod_boundary))
    times = [1e-6*k for k in range(1, int(end_us)+1)]
    r, refs = rod_reference(times)
    r2 = .5*(fluid.rf[:-1]**2+fluid.rf[1:]**2)
    rows, start = [], time.time()
    for t, a in zip(times, refs):
        model.advance(t)
        g = .5*(model.state.conserved[:, nz//2-1, 5]+model.state.conserved[:, nz//2, 5])
        # reference cell average of g = B_theta/r with weight r dr: integral of (a/r) dr over integral of r dr
        exact = np.array([np.trapezoid(np.interp(x, r, a)/x, x)/(.5*(hi*hi-lo*lo))
                          for lo, hi in zip(fluid.rf[:-1], fluid.rf[1:]) for x in [np.linspace(lo, hi, 201)]])
        d = model.diagnostics()
        work = d["toroidal_boundary_work_j"]
        rows.append(dict(time_s=model.time, max_error_over_k0=float(np.max(np.abs(r2*(g-exact)))/ROD_K0),
                         midplane_a_over_k0=[float(x) for x in (r2*g/ROD_K0)[::max(1, nr//6)]],
                         boundary_work_j=work, toroidal_magnetic_energy_j=d["toroidal_magnetic_energy_j"],
                         ohmic_heat_j=d["ohmic_heat_j"],
                         energy_theorem_residual_over_work=(d["toroidal_magnetic_energy_j"]+d["ohmic_heat_j"]
                                                            +d["toroidal_kinetic_energy_j"]-work)/work,
                         fluid_energy_residual_j=d["fluid_energy_balance_residual_j"],
                         response_energy_residual_j=d["response_energy_balance_residual_j"],
                         maximum_speed_m_s=d["maximum_speed_m_s"], steps=d["steps"]))
    return dict(case="rod", nr=nr, nz=nz, wall_s=time.time()-start, reference_intervals=4000, at=rows)


CASES = dict(torsional=torsional, rotor=rotor, alfvenic=alfvenic, whistler=whistler, rod=rod)

if __name__ == "__main__":
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    case = sys.argv[2]
    extra = dict(a.split("=", 1) for a in sys.argv[3:] if "=" in a)
    sizes = [int(a) for a in sys.argv[3:] if "=" not in a]
    results = []
    for nr in sizes:
        row = CASES[case](nr, **{k: float(v) for k, v in extra.items()})
        results.append(row)
        print(json.dumps(row), flush=True)
        (out/f"{case}{''.join(f'_{k}{v}' for k, v in extra.items())}.json").write_text(json.dumps(results, indent=1))
