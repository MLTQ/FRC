"""Verification of the 3D hybrid core (plasma_bench/hybrid3d.py) against independent exact results.

1. Boris orbit in a uniform B: gyration phase error against the exact orbit, which should be second order in
   Omega dt, and exact energy.
2. Parallel circularly polarized waves in a uniform cold plasma (cold ions as particles, massless cold
   electrons), against the exact cold-ion hybrid (Hall-MHD) dispersion
       omega^2 = k^2 v_A^2 (1 + omega/Omega_i)
   whose positive root is the whistler branch and whose negative root is the ion-cyclotron branch
   (|omega| -> Omega_i). The eigenmode is
       dB = A (cos phi, -sin phi, 0),  du = -Omega omega A/(k B0 (omega+Omega)) (cos phi, -sin phi, 0),
   with phi = k z - omega t. The frequency is measured from the phase of
   c(t) = <(Bx - i By) e^{-ikz}> = A e^{-i omega t}, at 16, 32 and 64 cells per wavelength.
   Also reported: div B and the drift of total energy (ion kinetic + magnetic) against the wave energy.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE, m_p as MP
from plasma_bench.hybrid3d import Box, Hybrid, HybridParams, boris

B0, N0, MASS = .1, 1e20, MP


def orbit_check():
    """Single ion gyrating in Bz = B0, E = 0: phase error after 10 periods against Omega dt."""
    omega = QE*B0/MASS
    rows = []
    for steps_per_period in (16, 32, 64, 128):
        dt = 2*np.pi/omega/steps_per_period
        x = np.array([[0., 0., 0.]]); v = np.array([[1e5, 0., 0.]])
        # start v at -dt/2: rotate the exact velocity back by half a step
        a = -.5*omega*dt
        v = np.array([[1e5*np.cos(a), -1e5*np.sin(a), 0.]])
        e = np.zeros((3, 4, 4, 4)); b = np.zeros((3, 4, 4, 4)); b[2] = B0
        lo, d, sh, ln = np.zeros(3), np.full(3, 1.), np.array([4, 4, 4]), np.full(3, 4.)
        x[:] = 2.
        n = 10*steps_per_period
        for _ in range(n):
            boris(x, v, e, b, QE/MASS, dt, lo, d, sh, ln)
        # velocity at t = n dt - dt/2: exact phase -omega (n dt - dt/2) (ions gyrate clockwise about +z)
        exact = -omega*(n*dt-.5*dt)
        phase = np.arctan2(v[0, 1], v[0, 0])
        err = (phase-exact+np.pi) % (2*np.pi)-np.pi
        rows.append(dict(steps_per_period=steps_per_period, phase_error_rad=float(err),
                         speed_relative_error=float(abs(np.linalg.norm(v)/1e5-1))))
    return rows


def wave(kdi, branch, cells, periods=2., ppc_z=8, amp=1e-3, cfl=.5):
    omega_i = QE*B0/MASS
    va = B0/np.sqrt(MU0*N0*MASS)
    di = va/omega_i
    k = kdi/di
    disc = np.sqrt(1+4*omega_i**2/(k*va)**2)
    omega = (k*va)**2/(2*omega_i)*(1+disc if branch == "whistler" else 1-disc)
    length = 2*np.pi/k
    nx = ny = 2
    dz = length/cells
    box = Box((nx, ny, cells), (0., 0., 0.), (nx*dz, ny*dz, length))
    a = amp*B0
    zf = (np.arange(cells)+.5)*dz
    bx = np.broadcast_to(a*np.cos(k*zf), (nx, ny, cells)).copy()
    by = np.broadcast_to(-a*np.sin(k*zf), (nx, ny, cells)).copy()
    bz = np.full((nx, ny, cells), B0)
    period = 2*np.pi/abs(omega)
    dt = min(.05/omega_i, period/64, .25*dz/(abs(omega)/k*(1+amp)+1e-30))
    # cold particles on a lattice: 1 x 1 x ppc_z per cell column
    zs = (np.arange(cells*ppc_z)+.5)*dz/ppc_z
    xs = (np.arange(nx)+.5)*dz; ys = (np.arange(ny)+.5)*dz
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    pos = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    phi = k*pos[:, 2]+omega*.5*dt
    du = -omega_i*omega*a/(k*B0*(omega+omega_i))
    vel = np.column_stack((du*np.cos(phi), -du*np.sin(phi), np.zeros(len(pos))))
    w = N0*box.cell_volume/(ppc_z)
    model = Hybrid(box, HybridParams(ion_mass=MASS, field_cfl=cfl), (bx, by, bz), pos, vel, w, dt, n_ref=N0)
    steps = int(np.ceil(periods*period/dt))
    times, amps, energy, divs = [], [], [], []
    kz = np.exp(-1j*k*zf)
    v_prev = model.v.copy()
    start = time.time()
    for s in range(steps+1):
        c = np.mean((model.b[0]-1j*model.b[1])*kz)
        times.append(model.time); amps.append(c)
        if s:
            kin = .5*(model.kinetic_energy(v_prev)+model.kinetic_energy())   # K at integer time ~ mean of half steps
            energy.append(kin+model.magnetic_energy())
        if s == steps:
            break
        v_prev = model.v.copy()
        model.step()
        divs.append(model.max_divergence())
    phase = np.unwrap(np.angle(amps))
    slope = np.polyfit(times, phase, 1)[0]
    measured = -slope
    vol = box.cell_volume*nx*ny*cells
    wave_energy = (a**2/(2*MU0))*vol*(1+(MU0*N0*MASS)*(du/a)**2)   # magnetic + kinetic of the eigenmode
    return dict(kdi=kdi, branch=branch, cells_per_wavelength=cells, omega_exact=omega, omega_measured=float(measured),
                relative_error=float(measured/omega-1), amplitude_ratio=float(abs(amps[-1])/abs(amps[0])),
                energy_drift_over_wave_energy=float((energy[-1]-energy[0])/wave_energy) if energy else None,
                max_div_b_relative=float(max(divs)), ion_steps=steps, field_substeps=int(model.substeps),
                dt_omega_i=float(dt*omega_i), wall_s=time.time()-start)


def main(out=Path("outputs/hybrid3d_verification")):
    out.mkdir(parents=True, exist_ok=True)
    result = dict(orbit=orbit_check())
    for r in result["orbit"]:
        print(json.dumps(r), flush=True)
    rows = []
    for kdi in (.25, .5, 1., 2.):
        for branch in ("whistler", "ion_cyclotron"):
            for cells in (16, 32, 64):
                rows.append(wave(kdi, branch, cells))
                print(json.dumps({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in rows[-1].items()}), flush=True)
    result["waves"] = rows
    (out/"results.json").write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main(*(sys.argv[1:2] and [Path(sys.argv[1])]))


def device_fast_mode(dx, dz=.05, periods=2., amp=.01, ppc=(2, 2, 2), n0=1e22, b0=.2, periodic_length=None, coil=.165):
    """Radial compressional (fast) eigenmode of a cold uniform plasma column in the 12" device: plasma for
    r < a (gas cells), exact vacuum for a < r < b (the gap), perfect conductor at b (the coil). With k_z = 0 the
    Hall term drops out (its E is a pure radial gradient) and the exact frequency is omega = x v_A/a with
        2 J1(x)/x + ((b^2 - a^2)/a^2) J0(x) = 0
    (flux conservation inside the conductor: the gap's uniform dBz carries the flux the plasma's J0 mode moves).
    a and b are the radii of equal staircase area at the midplane. The eigenmode is
        u_r = U J1(k r) cos(omega t),   dBz = -(B0 U k/omega) J0(k r) sin(omega t).
    Cold two-fluid theory gives omega = k v_A for k perpendicular to B at any omega/Omega_i, but the eigenmode then
    also carries an azimuthal ion velocity of relative size omega/Omega_i. The radial-only start above is the
    eigenmode only for omega << Omega_i, hence the high default density (omega/Omega_i ~ 0.05). At 1e20 m^-3
    (omega/Omega_i ~ 0.5) the same start also excites ion-cyclotron oscillations, and the fit is off by 5-9%."""
    from scipy.optimize import brentq, curve_fit
    from scipy.special import j0, j1
    from plasma_bench.hybrid3d_device import Device, device_box, cell_types, DeviceHybrid, GAS, CONDUCTOR, periodic
    base = Device(coil_radius=coil)
    if periodic_length:     # infinitely long column: the gap is a ring carrying the axial-flux harmonic
        dev, box = periodic(base), device_box(base, dx, dx, periodic_length=periodic_length)
    else:
        dev, box = base, device_box(base, dx, dz)
    types = cell_types(box, dev)
    mid = types[:, :, box.shape[2]//2]
    a = np.sqrt((mid == GAS).sum()*dx*dx/np.pi); b = np.sqrt((mid != CONDUCTOR).sum()*dx*dx/np.pi)
    g = (b**2-a**2)/a**2
    xr = brentq(lambda s: 2*j1(s)/s+g*j0(s), 2.0, 3.8317) if g > 0 else 3.8317059702075125
    k = xr/a
    va = b0/np.sqrt(MU0*n0*MASS)
    omega = k*va
    omega_i = QE*b0/MASS
    dt = min(.2/omega_i, 2*np.pi/omega/60)
    # cold ions on a lattice in every gas cell
    gi = np.argwhere(types == GAS)
    offs = np.array([((i+.5)/ppc[0], (j+.5)/ppc[1], (l+.5)/ppc[2]) for i in range(ppc[0]) for j in range(ppc[1]) for l in range(ppc[2])])
    pos = (np.asarray(box.lo)+(gi[:, None, :]+offs[None, :, :])*box.d).reshape(-1, 3)
    r = np.hypot(pos[:, 0], pos[:, 1])
    U = amp*va
    ur = U*j1(k*r)*np.cos(-.5*omega*dt)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx, sx = np.where(r > 0, pos[:, 0]/r, 1.), np.where(r > 0, pos[:, 1]/r, 0.)
    vel = np.column_stack((ur*cx, ur*sx, np.zeros(len(pos))))
    w = n0*box.cell_volume/len(offs)
    bf = [np.zeros(box.shape), np.zeros(box.shape), np.full(box.shape, b0)]
    model = DeviceHybrid(box, dev, HybridParams(ion_mass=MASS), bf, pos, vel, w, dt, n_ref=n0)
    # axis probe: the four Bz faces around the axis at the midplane
    cx0 = int(round(-box.lo[0]/box.d[0])); kz = slice(None) if periodic_length else box.shape[2]//2
    steps = int(np.ceil(periods*2*np.pi/omega/dt))
    t, sig, energy = [], [], []
    start = time.time()
    e0 = model.kinetic_energy()+model.magnetic_energy()
    for s in range(steps):
        model.step()
        t.append(model.time)
        sig.append(float(np.mean(model.b[2][cx0-1:cx0+1, cx0-1:cx0+1, kz])-b0))
        energy.append(model.kinetic_energy()+model.magnetic_energy())
    t, sig = np.array(t), np.array(sig)
    expected_amp = b0*U*k/omega
    fit, _ = curve_fit(lambda tt, A, om, ph: A*np.sin(om*tt+ph), t, sig, p0=(-expected_amp, omega, 0.))
    return dict(dx=dx, dz=dz, periodic=bool(periodic_length), coil_radius=coil, n0=n0, omega_over_omega_i=float(omega/omega_i), a_eff=float(a), b_eff=float(b), x_root=float(xr), omega_exact=float(omega),
                omega_measured=float(abs(fit[1])), relative_error=float(abs(fit[1])/omega-1),
                axis_amplitude_ratio=float(abs(fit[0])/expected_amp), steps=steps, field_substeps=int(model.substeps),
                particles_lost=float(model.lost["count"]), energy_change_relative=float((energy[-1]-e0)/e0),
                max_div_b=model.max_divergence(), wall_s=time.time()-start)


def rigid_rotor(dx, length=.2, ppc=(2, 2, 2), t_end=3e-6, b_e=.2, r_s=.11, kappa=1., ti_ev=75., te_ev=75.,
                floor_frac=.02, seed=0):
    """Exact kinetic equilibrium of the hybrid model in an infinitely long (z-periodic) 12" tube: non-rotating
    Maxwellian ions (T_i) and isothermal massless electrons (T_e) rotating rigidly at omega_e, with
        u = kappa (2 r^2/r_s^2 - 1),  B_z = B_e tanh u,  n = n_m sech^2 u,  n_m (T_i + T_e) = B_e^2/(2 mu0),
        omega_e = 4 B_e kappa / (mu0 e n_m r_s^2),  psi = B_e r_s^2/(4 kappa) ln(cosh u / cosh kappa).
    Ions: f(H) with H = m v^2/2 + e phi is stationary in any static field, and the electron force balance closes
    with ln n = e omega_e psi/(T_i + T_e) + const (a Liouville equation in r^2). Nothing should evolve except
    particle noise, the discretization of the edge, the declared density floor in Ohm's law where n < floor,
    and absorption of the few ions at the wall. B comes from A_theta = psi/r on edges (exactly div-free)."""
    from plasma_bench.hybrid3d import curl_edges
    from plasma_bench.hybrid3d_device import Device, device_box, cell_types, DeviceHybrid, GAS, periodic
    dev = periodic(Device())
    box = device_box(Device(), dx, dx, periodic_length=length)
    types = cell_types(box, dev)
    n_m = b_e**2/(2*MU0*(ti_ev+te_ev)*QE)
    def u(r): return kappa*(2*r**2/r_s**2-1)
    def psi(r): return b_e*r_s**2/(4*kappa)*np.log(np.cosh(u(r))/np.cosh(kappa))
    lo, d = np.asarray(box.lo), box.d
    xs = [lo[a]+d[a]*np.arange(box.shape[a]) for a in range(3)]
    xh = [lo[a]+d[a]*(np.arange(box.shape[a])+.5) for a in range(3)]
    def a_theta_over_r(x, y):
        r2 = x**2+y**2
        rr = np.sqrt(r2)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(r2 > 0, psi(rr)/r2, b_e*np.tanh(-kappa)/2)    # psi/r^2 -> B_z(0)/2 on the axis
    X, Y, _ = np.meshgrid(xh[0], xs[1], xs[2], indexing="ij")             # Ax at (i+1/2, j, k)
    ax = -Y*a_theta_over_r(X, Y)
    X, Y, _ = np.meshgrid(xs[0], xh[1], xs[2], indexing="ij")             # Ay at (i, j+1/2, k)
    ay = X*a_theta_over_r(X, Y)
    bf = list(curl_edges(ax, ay, np.zeros(box.shape), d))
    rng = np.random.default_rng(seed)
    gi = np.argwhere(types == GAS)
    offs = np.array([((i+.5)/ppc[0], (j+.5)/ppc[1], (l+.5)/ppc[2]) for i in range(ppc[0]) for j in range(ppc[1]) for l in range(ppc[2])])
    pos = (lo+(gi[:, None, :]+offs[None, :, :])*d).reshape(-1, 3)
    r = np.hypot(pos[:, 0], pos[:, 1])
    w = n_m/np.cosh(u(r))**2*box.cell_volume/len(offs)
    vel = rng.normal(scale=np.sqrt(ti_ev*QE/MASS), size=pos.shape)
    omega_i = QE*b_e/MASS
    dt = .1/omega_i
    params = HybridParams(ion_mass=MASS, te_ev=te_ev, n_floor=floor_frac*n_m)
    model = DeviceHybrid(box, dev, params, bf, pos, vel, w, dt, n_ref=n_m)
    edges = np.arange(0, dev.wall_radius+d[0], d[0])
    def profiles():
        zc = xh[2]
        Xf, Yf = np.meshgrid(xh[0], xh[1], indexing="ij")
        rf = np.hypot(Xf, Yf)
        bz = model.b[2].mean(axis=2)
        idx = np.digitize(rf, edges)-1
        nb = len(edges)-1
        prof_b = np.array([bz[idx == q].mean() if (idx == q).any() else np.nan for q in range(nb)])
        Xn, Yn = np.meshgrid(xs[0], xs[1], indexing="ij")
        rn = np.hypot(Xn, Yn)
        nn = model.n.mean(axis=2)
        fx, fy = model.flux_half[0].mean(axis=2), model.flux_half[1].mean(axis=2)
        with np.errstate(invalid="ignore", divide="ignore"):
            uth = np.where(nn > 0, (-Yn*fx+Xn*fy)/np.maximum(rn, 1e-12)/nn, 0.)
        idn = np.digitize(rn, edges)-1
        prof_n = np.array([nn[idn == q].mean() if (idn == q).any() else np.nan for q in range(nb)])
        prof_u = np.array([uth[idn == q].mean() if (idn == q).any() else np.nan for q in range(nb)])
        return prof_b, prof_n, prof_u
    centres = .5*(edges[:-1]+edges[1:])
    b0, n0, u0 = profiles()
    e0 = model.kinetic_energy()+model.magnetic_energy()
    start = time.time()
    steps = int(np.ceil(t_end/dt))
    for _ in range(steps):
        model.step()
    b1, n1, u1 = profiles()
    inside = centres < .9*dev.wall_radius
    exact_b = b_e*np.tanh(u(centres)); exact_n = n_m/np.cosh(u(centres))**2
    return dict(dx=dx, ppc=int(len(offs)), particles=int(len(pos)), steps=steps, t_end_s=float(model.time),
                substeps_per_step=float(model.substeps/steps), floor_over_nm=floor_frac,
                initial_bz_error_vs_exact=float(np.nanmax(np.abs(b0-exact_b)[inside])/b_e),
                bz_change_over_be=float(np.nanmax(np.abs(b1-b0)[inside])/b_e),
                n_change_over_nm=float(np.nanmax(np.abs(n1-n0)[inside])/n_m),
                initial_n_error_vs_exact=float(np.nanmax(np.abs(n0-exact_n)[inside])/n_m),
                axis_bz_start_end=[float(b0[0]), float(b1[0])], max_ion_rotation_over_vti=float(np.nanmax(np.abs(u1[inside]))/np.sqrt(ti_ev*QE/MASS)),
                lost_fraction=float(model.lost["count"]/w.sum()), energy_change_relative=float((model.kinetic_energy()+model.magnetic_energy()-e0)/e0),
                max_div_b=model.max_divergence(), wall_s=time.time()-start,
                profiles=dict(r=centres.tolist(), bz_start=b0.tolist(), bz_end=b1.tolist(), n_start=n0.tolist(), n_end=n1.tolist(),
                              u_theta_end=u1.tolist(), bz_exact=exact_b.tolist(), n_exact=exact_n.tolist()))
