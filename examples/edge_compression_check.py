"""Adiabatic theta-pinch compression of a magnetized column: an analytic plasma-vacuum edge check (frc-ibl.11.9).

    python -m examples.edge_compression_check OUT model=floor|mask|edge|ledger [ppc=8 seed=0 t_ramp_us=4 t_end_us=5 min_ions=1 max_scale=2 electron_edge=1]

A uniform hydrogen column (radius a0, density n0, Ti = Te, internal field B_in0) sits in the 12" device inside an
exact-vacuum annulus, in pressure balance with the external field B_out0. A uniform solenoid field B_s(t) then
rises smoothly to B_s1 and holds. The column ends are the tube's absorbing end walls and the end planes are
conductors, so the flux through every cross-section inside the coil is fixed by the end-plane data:
    B_out = (B_s b^2 - B_out0 a0^2)/(b^2 - a^2)          (b = coil radius of equal staircase area)
Quasi-static, ideal, sharp-boundary compression conserves the column flux (B_in a^2), the ion magnetic moment
(collisionless double-adiabatic, T_perp ~ B, T_par fixed: p_perp a^4 fixed) and the electron adiabat
(gamma = 5/3, p_e a^(10/3) fixed). The radius follows from radial pressure balance
    p_i0 x^4 + p_e0 x^(10/3) + (B_in0 x^2)^2/2mu0 = B_out^2/2mu0,     x = a0/a,
and no energy is dissipated. The measurement compares the central ion radius sqrt(2<r^2>) (|z| < 0.3 m) and the
core Bz (r < 6 cm, |z| < 0.3 m)
(both relative to their own t = 0 values) with this curve, and books every energy term, including all numerical
removals, against the measured source work. Models:
- floor: DeviceHybrid; the empty gas annulus is the 1 eV, floor-density Ohm's-law medium;
- mask:  VacuumMaskHybrid (cells with no macro-ion are exact vacuum);
- edge:  PlasmaEdgeHybrid (mask plus an interface floor taken on the plasma side; with electron_edge=1 also the
         EdgeEnergyElectrons closure edge);
- ledger: CarrierLedgerHybrid (the floor field closure with a carrier-class heat ledger; electron_edge=1 adds
          the electron-edge transport without any field mask).
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE, m_p
from scipy.optimize import brentq
from plasma_bench.hybrid3d import HybridParams, curl_edges, curl_faces
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, DeviceHybrid, GAS, CONDUCTOR
from plasma_bench.hybrid3d_vacuum_mask import VacuumMaskHybrid
from plasma_bench.hybrid3d_edge import PlasmaEdgeHybrid
from plasma_bench.hybrid3d_carrier_ledger import CarrierLedgerHybrid

MASS = m_p+9.1093837e-31
CASE = dict(a0=.12, n0=5e19, ti_ev=20., te_ev=20., b_out0=.05, b_s1=.15, floor=.05, dx=.02, dz=.04, dt_ns=2.)


class UniformSolenoid:
    """Prescribed uniform Bz = B_s(t) inside the coil: edge A = B_s (-y/2, x/2, 0), faces B = curl A (exact)."""
    def __init__(self, box, waveform, derivative):
        xs = [box.lo[a]+box.d[a]*np.arange(box.shape[a]) for a in range(3)]
        xh = [x+.5*box.d[a] for a, x in enumerate(xs)]
        ax = -.5*np.meshgrid(xh[0], xs[1], xs[2], indexing='ij')[1]
        ay = .5*np.meshgrid(xs[0], xh[1], xs[2], indexing='ij')[0]
        self.edges = np.array([[ax, ay, np.zeros(box.shape)]])
        self.faces = np.array([curl_edges(*self.edges[0], box.d)])
        self.current_basis = np.array([np.asarray(curl_faces(*self.faces[0], box.d))/MU0])
        self.waveform, self.derivative = waveform, derivative

    def coefficients(self, time_s):
        return np.array([self.waveform(time_s)])

    def __call__(self, time_s):
        return self.waveform(time_s)*self.faces[0]

    def edge_derivative(self, time_s):
        return self.derivative(time_s)*self.edges[0]

    def edge_current(self, time_s):
        return self.waveform(time_s)*self.current_basis[0]


def analytic(case, b_s, b_coil):
    """Quasi-static radius ratio a/a0 and B_in/B_in0 for source field b_s."""
    a0, b0 = case['a0'], case['b_out0']
    pi0, pe0 = case['n0']*case['ti_ev']*QE, case['n0']*case['te_ev']*QE
    bin0 = np.sqrt(b0**2-2*MU0*(pi0+pe0))

    def balance(x):
        a = a0/x
        bout = (b_s*b_coil**2-b0*a0**2)/(b_coil**2-a**2)
        return pi0*x**4+pe0*x**(10/3)+(bin0*x**2)**2/(2*MU0)-bout**2/(2*MU0)
    x = 1. if balance(1.) >= 0 else brentq(balance, 1., 20.)
    return 1/x, x**2, bin0


def run(out, model='edge', ppc=8, seed=0, t_ramp_us=4., t_end_us=5., min_ions=1, max_scale=2., electron_edge=1, every_ns=50.):
    case = dict(CASE)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    dev = Device()
    box = device_box(dev, case['dx'], case['dz'])
    types = cell_types(box, dev)
    mid = types[:, :, box.shape[2]//2]
    b_coil = float(np.sqrt((mid != CONDUCTOR).sum()*case['dx']**2/np.pi))
    t_ramp, b0, b1 = t_ramp_us*1e-6, case['b_out0'], case['b_s1']

    def waveform(t):
        return b0+(b1-b0)*.5*(1-np.cos(np.pi*min(max(t, 0.), t_ramp)/t_ramp))

    def derivative(t):
        return (b1-b0)*.5*np.pi/t_ramp*np.sin(np.pi*t/t_ramp) if 0 <= t <= t_ramp else 0.
    source = UniformSolenoid(box, waveform, derivative)
    _, _, bin0 = analytic(case, b0, b_coil)
    # ions: a cubic lattice in gas cells, kept where r < a0; Maxwellian, no drift
    side = int(round(ppc**(1/3)))
    gi = np.argwhere(types == GAS)
    offs = np.array([((i+.5)/side, (j+.5)/side, (l+.5)/side) for i in range(side) for j in range(side) for l in range(side)])
    pos = (np.asarray(box.lo)+(gi[:, None, :]+offs[None, :, :])*box.d).reshape(-1, 3)
    pos = pos[np.hypot(pos[:, 0], pos[:, 1]) < case['a0']]
    rng = np.random.default_rng(int(seed))
    vel = rng.normal(size=pos.shape)*np.sqrt(case['ti_ev']*QE/MASS)
    w = case['n0']*box.cell_volume/side**3
    # field: source plus a straight diamagnetic tube B_in0 - B_out0 inside a0 (exactly solenoidal)
    xs = [box.lo[a]+box.d[a]*(np.arange(box.shape[a])+.5) for a in range(2)]
    inside = (np.hypot(xs[0][:, None], xs[1][None, :]) < case['a0'])[:, :, None]*np.ones(box.shape)
    bf = [c.copy() for c in source(0.)]
    bf[2] = bf[2]+(bin0-b0)*inside
    params = HybridParams(ion_mass=MASS, te_ev=case['te_ev'], gamma_e=5/3, n_floor=case['floor']*case['n0'])
    thermal = dict(closure='energy', gamma=5/3, resistivity='spitzer', coulomb_log=10., transport_cfl=.8)
    cls = dict(floor=DeviceHybrid, mask=VacuumMaskHybrid, edge=PlasmaEdgeHybrid, ledger=CarrierLedgerHybrid)[model]
    extra = (dict(electron_edge=bool(electron_edge)) if model == 'ledger' else {} if model == 'floor'
             else dict(min_ions=int(min_ions)))
    if model == 'edge':
        extra.update(max_scale=max_scale, electron_edge=bool(electron_edge))
    m = cls(box, dev, params, bf, pos, vel, w, case['dt_ns']*1e-9, n_ref=case['n0'], source=source, thermal=thermal, **extra)
    zc = box.lo[2]+box.d[2]*np.arange(box.shape[2])          # Bz faces sit on z nodes
    core = ((np.hypot(xs[0][:, None], xs[1][None, :]) < .06)[:, :, None]) & (np.abs(zc) < .3)[None, None, :]

    def linkage():
        jp = np.asarray(curl_faces(*(np.asarray(m.b)-source(m.time)), m.d))/MU0*np.asarray(m.edge_mask)
        return float(np.einsum('axyz,axyz->', jp, source.edges[0])*box.cell_volume)

    def radius():
        sel = np.abs(m.x[:, 2]) < .3
        return float(np.sqrt(2*np.mean(m.x[sel, 0]**2+m.x[sel, 1]**2)))

    def axis_field():
        return float(np.mean(m.b[2][core]))          # core Bz: r < 6 cm, |z| < 30 cm

    def accounted():
        e = m.electrons
        return (m.response_magnetic_energy()+m.kinetic_energy()+e.energy(m.n)+m.lost['energy_j']+e.absorbed_energy_j
                +e.cutoff_removed_heat_j+getattr(m, 'mask_projection_response_j', 0.))
    a_start, b_start, e_start = radius(), axis_field(), accounted()
    work, phi, b_prev = 0., linkage(), waveform(0.)
    rows, start, next_row = [], time.time(), 0.

    def record():
        ratio, bratio, _ = analytic(case, waveform(m.time), b_coil)
        e = m.electrons
        row = dict(time_us=m.time*1e6, b_source_t=waveform(m.time), radius_ratio=radius()/a_start,
                   radius_ratio_analytic=ratio, axis_b_ratio=axis_field()/b_start, axis_b_ratio_analytic=bratio,
                   response_magnetic_j=m.response_magnetic_energy(), kinetic_j=m.kinetic_energy(),
                   electron_internal_j=e.energy(m.n), lost_ion_energy_j=m.lost['energy_j'],
                   ohmic_heat_j=e.ohmic_heat_j, cutoff_removed_heat_j=e.cutoff_removed_heat_j,
                   empty_node_energy_j=e.vacuum_removed_energy_j, boundary_enthalpy_j=e.enthalpy_outflow_j,
                   source_work_j=work, energy_residual_j=accounted()-e_start-work, substeps=m.substeps,
                   electron_pressure_centering_residual_j=e.pressure_centering_residual_j,
                   electron_grid_residual_j=e.energy_balance_residual_j, electron_roundoff_j=e.roundoff_energy_j,
                   wall_s=time.time()-start)
        if hasattr(m, 'mask_diagnostics'):
            row.update(m.mask_diagnostics())
        if hasattr(m, 'ledger_diagnostics'):
            row.update(m.ledger_diagnostics())
        rows.append(row)
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    while m.time < t_end_us*1e-6-1e-15:
        if m.time >= next_row-1e-15:
            record()
            next_row += every_ns*1e-9
        m.step()
        phi_new, b_new = linkage(), waveform(m.time)
        work -= .5*(phi+phi_new)*(b_new-b_prev)
        phi, b_prev = phi_new, b_new
    record()
    info = dict(case=case, model=model, ppc=side**3, seed=seed, min_ions=min_ions, max_scale=max_scale if model == 'edge' else None,
                electron_closure=type(m.electrons).__name__, t_ramp_us=t_ramp_us, t_end_us=t_end_us,
                coil_radius_equal_area_m=b_coil, b_in0_t=float(bin0), particles=int(len(pos)), box=list(box.shape),
                initial_radius_m=a_start)
    (out/'edge_compression.json').write_text(json.dumps(dict(info=info, rows=rows), indent=1))
    return rows


if __name__ == '__main__':
    args = dict(a.split('=', 1) for a in sys.argv[2:])
    run(sys.argv[1], **{k: (v if k == 'model' else float(v)) for k, v in args.items()})
