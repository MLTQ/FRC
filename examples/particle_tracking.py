"""3D test-particle orbits in saved FRC formation runs (plasma_bench/particle_tracking.py).

    python -m examples.particle_tracking verify OUT
    python -m examples.particle_tracking ions RUN OUT [start_us=10] [end_us=30] [count=4000] [seed=1]
    python -m examples.particle_tracking rmf RUN OUT [frequency=1e6] [amplitude=5e-3] [parity=even] [count=300]
    python -m examples.particle_tracking multipoles RUN OUT [cases=none,q3,q5,q7,q9,o7] [count=4000] [seed=1]

verify: analytic checks of the pusher and fields (gyration, ExB and grad-B drifts, mirror loss cone, coil
fields) and energy / p_theta conservation with convergence order in the saved FRC fields of a run
(RUN defaults to outputs/fine_uniform_h2).
ions: a thermal ion population loaded inside the separatrix of RUN at start_us, tracked to end_us without
and with the straight quadrupole at 7% of B_e (set at the wall and at the separatrix): orbit classes,
confinement against time, end/wall loss maps, s-bar. Writes OUT/summary.json, figures, and
RUN/particles.json (decimated orbits for examples/formation_display.py).
rmf: electrons in the same fields plus an RMF antenna, with a no-RMF control: azimuthal entrainment.
multipoles: the same ions with static multipoles of different order and strength (q7 = 4-bar quadrupole at 7%
of B_e, o7 = 8-bar octupole at 7%, mean |B| on the wall circle): confined fraction, wall/end split and the
m = 2 / m = 4 loss-sector amplitude, paired against the no-multipole case. Cases run in parallel processes.

Fields (companion particle_tracking.md): the saved gas-node flux, extended into the vacuum exactly as in
examples/formation_display.py (harmonic extension of the induced flux plus the section coil flux of the
saved currents), splined as psi = r^2 g, linear in time between the 0.1 us saved states;
E_theta = -(dpsi/dt)/r. The single-fluid model has no toroidal field or flow, so its -v x B is purely
azimuthal and already contained in the Faraday E_theta; E_r = E_z = 0. velocity.npz is used for the
loaded ions' mean flow and for a consistency check of E_theta against -v x B.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.constants import e as QE, mu_0 as MU0, m_e as ME
from scipy.special import ellipe, ellipk
from plasma_bench.particle_tracking import (SPECIES, SURFACES, FluxField, Superposition, UniformField, GradientField,
                                            Walls, Population, track, canonical_angular_momentum, maxwellian,
                                            momentum_from_velocity, quadrupole, multipole, rmf_antenna, gamma_of)
from examples.formation_display import case_from_setup, geometry, fields as grid_fields
from examples.nucte_formation import coil_sections, ION_MASS
from examples.nucte_compare import load, closed_region, annular_closed_region, plasma_inside

ION = {"H2": "H+", "D2": "D+"}
CLASSES = ("betatron", "mixed", "drift", "cyclotron", "lost in window")
COLORS = {"betatron": "#ffffff", "mixed": "#d6b8ff", "drift": "#ff7ad9", "cyclotron": "#8fe3ff", "lost in window": "#ff9f43"}
QUAD_FRACTION = 0.07


# ---------------------------------------------------------------- saved run

class Run:
    """Saved formation run: geometry, frames and the FluxField over a time window."""

    def __init__(self, run, start_s, end_s):
        self.path = Path(run)
        self.setup, self.rf, self.zf, frames, self.history = load(run)
        c = self.case = case_from_setup(self.setup)
        self.i0, self.r_in, fluid, self.grid, mapping = geometry(c)
        self.walls = Walls(c.wall_radius, self.r_in, c.half_length)
        self.ion = ION[c.gas]
        self.density_mass = ION_MASS[c.gas]
        times = np.array([self.history[k]["time_s"] for k, *_ in frames])
        pick = [n for n, t in enumerate(times) if start_s-1e-12 <= t <= end_s+1e-12]
        if len(pick) < 2:
            raise ValueError(f"{run}: fewer than two saved states in [{start_s}, {end_s}] s")
        self.frames = [frames[n] for n in pick]
        self.times = times[pick]
        z0 = mapping.z_offset
        rows, cols = slice(self.i0, self.i0+c.nr+1), slice(z0, z0+c.nz+1)
        coil = coil_sections(c)
        rr, zz = np.meshgrid(self.grid.r, self.grid.z, indexing="ij")
        unit = np.array([f[0] for f in coil.field(rr, zz)])
        psis, self.reproduction = [], 0.
        for k, psi_gas, rho, p in self.frames:
            source = np.tensordot(np.asarray(self.history[k]["circuit"]["vacuum_current_a"]), unit, axes=1)
            induced = np.zeros(self.grid.shape)
            induced[rows, cols] = psi_gas-source[rows, cols]
            psi = mapping.harmonic_extension(induced)+source
            self.reproduction = max(self.reproduction, float(abs(psi[rows, cols]-psi_gas).max()/abs(psi_gas).max()))
            psis.append(psi)
        self.psi_full = np.array(psis)
        self.field = FluxField(self.grid.r, self.grid.z, self.times, self.psi_full)
        velocity = np.load(self.path/"velocity.npz")
        self.velocity = {k: velocity[f"t{k:04d}"].reshape(2, c.nr, c.nz) for k, *_ in self.frames}

    def region(self, n):
        k, psi, rho, p = self.frames[n]
        region = annular_closed_region(self.rf, self.zf, psi) if self.rf[0] > 0 else closed_region(self.rf, self.zf, psi)
        if region is None:
            return None
        nbar, temp = plasma_inside(self.rf, self.zf, rho, p, region["profile"], region["lo"], region["hi"],
                                   self.density_mass, region.get("inner_profile"))
        region.update(density_m3=nbar, total_temperature_ev=temp,
                      separatrix_psi=region.get("separatrix_flux_wb", 0.)/(2*np.pi))
        return region

    def index(self, t):
        return int(np.argmin(abs(self.times-t)))

    def external_field(self, t):
        """|B_z| at the wall, midplane: the field outside the FRC (B_e)."""
        X = np.array([[self.case.wall_radius*(1-1e-9), 0., 0.]])
        return float(abs(self.field.fields(X, t)[0][0, 2]))

    def load(self, count, species, t, rng):
        """Thermal particles inside the separatrix at saved state t: cells chosen with probability n r dr dz,
        positions uniform in volume within a cell and kept only where psi < psi_separatrix; local Maxwellian at
        T = p/(2 n) (Te = Ti) shifted by the cell's fluid velocity (v_r, v_z)."""
        n_state = self.index(t)
        k, psi, rho, p = self.frames[n_state]
        region = self.region(n_state)
        if region is None:
            raise ValueError(f"no closed flux at {t*1e6:.1f} us")
        rf, zf = self.rf, self.zf
        rc = .5*(rf[:-1]+rf[1:])
        psi_c = .25*(psi[:-1, :-1]+psi[1:, :-1]+psi[:-1, 1:]+psi[1:, 1:])
        inside = psi_c < region["separatrix_psi"]
        inside[:, :region["lo"]] = False
        inside[:, region["hi"]:] = False
        density = rho/self.density_mass
        weight = density*rc[:, None]*inside
        temperature = p/(2*density*QE)
        vr, vz = self.velocity[k]
        mass = SPECIES[species][0]
        X, V = [], []
        while sum(len(x) for x in X) < count:
            cells = rng.choice(weight.size, size=2*count, p=(weight/weight.sum()).ravel())
            i, j = np.unravel_index(cells, weight.shape)
            r = np.sqrt(rf[i]**2+rng.random(cells.size)*(rf[i+1]**2-rf[i]**2))
            z = zf[j]+rng.random(cells.size)*(zf[j+1]-zf[j])
            phi = 2*np.pi*rng.random(cells.size)
            x = np.column_stack((r*np.cos(phi), r*np.sin(phi), z))
            ok = self.field.psi(x, self.times[n_state]) < region["separatrix_psi"]
            v = maxwellian(rng, cells.size, temperature[i, j], mass)
            v[:, 0] += vr[i, j]*np.cos(phi)
            v[:, 1] += vr[i, j]*np.sin(phi)
            v[:, 2] += vz[i, j]
            X.append(x[ok]), V.append(v[ok])
        X, V = np.concatenate(X)[:count], np.concatenate(V)[:count]
        return Population(species, X, momentum_from_velocity(V)), region


def multipole_per_ampere(radius, order=2, reduce=np.mean):
    """Mean (or `reduce`) |B| [T/A] of the 2*order-bar multipole (order 2: the dossier quadrupole) around the
    midplane circle of `radius`. The 7% B_e level is set on this mean; the quadrupole field varies x2 around
    the azimuth at the wall."""
    th = np.linspace(0, 2*np.pi, 18*order+1)[:-1]
    X = np.column_stack((radius*np.cos(th), radius*np.sin(th), 0*th))
    return float(reduce(np.linalg.norm(multipole(1., order).fields(X, 0.)[0], axis=1)))


def time_step(fields, walls, t_list, species, omega_dt, snapshot_s):
    """dt with omega_c,max dt <= omega_dt (|B| sampled on chamber points at the given times) and an integer
    number of steps per saved-state interval, so no kick straddles a change of dpsi/dt."""
    r = np.linspace(max(walls.r_inner, 1e-4), walls.r_outer*(1-1e-6), 40)
    z = np.linspace(-walls.z_end, walls.z_end, 161)
    phi = np.arange(8)*np.pi/4
    R, Z, P = np.meshgrid(r, z, phi, indexing="ij")
    X = np.column_stack(((R*np.cos(P)).ravel(), (R*np.sin(P)).ravel(), Z.ravel()))
    bmax = max(float(np.linalg.norm(fields.fields(X, t)[0], axis=1).max()) for t in t_list)
    mass, charge = SPECIES[species]
    target = omega_dt*mass/(abs(charge)*bmax)
    n = int(np.ceil(snapshot_s/target))
    return snapshot_s/n, bmax


# ---------------------------------------------------------------- verification

def circle_radius(xy):
    """Algebraic least-squares circle fit (exact for points on a circle): x^2 + y^2 = 2 a x + 2 b y + c."""
    A = np.column_stack((2*xy[:, 0], 2*xy[:, 1], np.ones(len(xy))))
    a, b, c = np.linalg.lstsq(A, (xy**2).sum(axis=1), rcond=None)[0]
    return float(np.sqrt(c+a*a+b*b))


def drift_fit(t, y, phase):
    """Guiding-centre drift: least squares y = a + v t + c cos(phase) + d sin(phase), with the measured gyro
    phase of the momentum, so the gyration is removed exactly for a constant gyro-frequency."""
    A = np.column_stack((np.ones_like(t), t, np.cos(phase), np.sin(phase)))
    return float(np.linalg.lstsq(A, y, rcond=None)[0][1])


def gyration_check():
    """Uniform B = 0.2 T: gyro-frequency (relativistic, omega = qB/(gamma m)) and radius u/(gamma omega) against
    theory; the Boris rotation angle 2 atan(omega dt/2) per step is the expected frequency error."""
    rows = []
    for species, energy_ev in (("H+", 100.), ("e-", 100.), ("e-", 5e5)):
        mass, charge = SPECIES[species]
        b0 = 0.2
        gamma = 1+energy_ev*QE/(mass*2.99792458e8**2)
        u = 2.99792458e8*np.sqrt(gamma**2-1)
        omega = abs(charge)*b0/(gamma*mass)
        radius = u/(gamma*omega)
        for omega_dt in (0.2, 0.1, 0.05):
            dt = omega_dt/omega
            steps = int(round(20*2*np.pi/omega_dt))
            pop = Population(species, np.zeros((1, 3)), np.array([[u, 0., 0.]]))
            e0 = pop.energy_ev()[0]
            phase, xs = [], []
            track(pop, UniformField((0, 0, b0)), 0., dt, steps,
                  observe=lambda t, p, B: (phase.append(np.arctan2(p.u[0, 1], p.u[0, 0])), xs.append(p.x[0, :2].copy())))
            measured = abs(np.unwrap(np.r_[0., phase])[-1])/(steps*dt)
            boris = 2*np.arctan(omega_dt/2)/dt
            rows.append(dict(species=species, energy_ev=energy_ev, gamma=gamma, omega_dt=omega_dt,
                             frequency_error=measured/omega-1, boris_prediction=boris/omega-1,
                             radius_error=circle_radius(np.array(xs))/radius-1,
                             energy_relative_change=float(pop.energy_ev()[0]/e0-1)))
    return rows


def drift_check():
    """E x B (uniform, E = 10 kV/m x, B = 0.2 T z) and grad-B (B = B0 (1 + x/L) z, L = 1 m) drift velocities
    from a fit that removes the gyration, against E x B/B^2 and (v_perp^2/(2 omega_c L)) (B0/B_gc)^2 with B at
    the guiding centre (residual O((rho/L)^2))."""
    rows = []
    b0, e0 = 0.2, 1e4
    for species in ("H+", "e-"):
        mass, charge = SPECIES[species]
        omega = abs(charge)*b0/mass
        for omega_dt in (0.5, 0.1):
            dt = omega_dt/omega
            steps = int(round(50*2*np.pi/omega_dt))
            pop = Population(species, np.zeros((1, 3)), np.array([[3e4, 0., 0.]]))
            ys, ph = [], []
            # gyration about the drift: phase of u - v_E
            track(pop, UniformField((0, 0, b0), (e0, 0, 0)), 0., dt, steps,
                  observe=lambda t, p, B: (ys.append(p.x[0, 1]), ph.append(np.arctan2(p.u[0, 1]+e0/b0, p.u[0, 0]))))
            t = dt*np.arange(1, steps+1)
            v = drift_fit(t, np.array(ys), np.unwrap(ph))
            rows.append(dict(test="ExB", species=species, omega_dt=omega_dt, measured_vy=v, analytic_vy=-e0/b0,
                             relative_error=v/(-e0/b0)-1))
    mass, charge = SPECIES["H+"]
    length, vperp = 1., 1e5
    omega = charge*b0/mass
    # launched at x = 0 moving +y: the guiding centre sits at x = rho, where B = B0 (1 + rho/L)
    analytic = vperp**2/(2*omega*length)/(1+vperp/omega/length)**2
    for omega_dt in (0.4, 0.2, 0.1, 0.05):
        dt = omega_dt/omega
        steps = int(round(200*2*np.pi/omega_dt))
        pop = Population("H+", np.zeros((1, 3)), np.array([[0., vperp, 0.]]))
        ys, ph = [], []
        track(pop, GradientField(b0, length), 0., dt, steps,
              observe=lambda t, p, B: (ys.append(p.x[0, 1]), ph.append(np.arctan2(p.u[0, 1], p.u[0, 0]))))
        t = dt*np.arange(1, steps+1)
        measured = drift_fit(t, np.array(ys), np.unwrap(ph))
        rows.append(dict(test="gradB", species="H+", omega_dt=omega_dt, measured_vy=measured, analytic_vy=analytic,
                         relative_error=measured/analytic-1, rho_over_L=vperp/omega/length))
    return rows


def loop_flux(a, current, r, z):
    """Exact flux per radian r A_phi of a circular loop of radius a at z = 0 (elliptic integrals)."""
    r = np.maximum(r, 1e-12)
    k2 = 4*a*r/((a+r)**2+z**2)
    k = np.sqrt(k2)
    a_phi = MU0*current/(np.pi*k)*np.sqrt(a/r)*((1-k2/2)*ellipk(k2)-ellipe(k2))
    return r*a_phi


def mirror_check(rng, count=2000):
    """Two coaxial loops (radius a = 0.2 m at z = +/-0.3 m): their exact flux on a grid -> FluxField (the
    production field path). Isotropic 10 eV protons launched within 5 mm of the axis at the midplane
    (rho/L ~ 0.015) are lost iff they pass the field maximum. Adiabatic prediction: lost iff
    sin^2(pitch) < B(x0)/B_max; isotropic loss fraction 1 - sqrt(1 - 1/R)."""
    a, L, current = 0.2, 0.3, 9.3e4
    r = np.linspace(0, 0.15, 61)
    z = np.linspace(-0.5, 0.5, 401)
    R, Z = np.meshgrid(r, z, indexing="ij")
    field = FluxField(r, z, [0.], loop_flux(a, current, R, Z-L)+loop_flux(a, current, R, Z+L))
    zz = np.linspace(-0.5, 0.5, 4001)
    b_axis = MU0*current*a**2/2*((a**2+(zz-L)**2)**-1.5+(a**2+(zz+L)**2)**-1.5)
    z_peak, bmax = float(abs(zz[np.argmax(b_axis)])), float(b_axis.max())
    mass, charge = SPECIES["H+"]
    speed = np.sqrt(2*10.*QE/mass)
    d = rng.normal(size=(count, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    rr = 0.005*np.sqrt(rng.random(count)); ph = 2*np.pi*rng.random(count)
    X = np.column_stack((rr*np.cos(ph), rr*np.sin(ph), np.zeros(count)))
    B0 = field.fields(X, 0.)[0]
    b0 = np.linalg.norm(B0, axis=1)
    b_centre = float(MU0*current*a**2/(a**2+L**2)**1.5)
    sin2 = 1-(np.einsum("nk,nk->n", d, B0)/b0)**2
    predicted = sin2 < b0/bmax
    pop = Population("H+", X, speed*d)
    dt = 0.1/(charge*bmax/mass)
    bounce = 4*z_peak/speed
    track(pop, field, 0., dt, int(1.5*bounce/dt), walls=Walls(0.15, 0., z_peak+0.1))
    lost = np.isfinite(pop.lost_time)
    ratio = bmax/b_centre
    wrong = lost != predicted
    margin = np.degrees(np.abs(np.arcsin(np.sqrt(sin2))-np.arcsin(np.sqrt(b0/bmax))))
    return dict(mirror_ratio=ratio, centre_field_t=b_centre, centre_field_spline_t=float(b0.mean()),
                rho_over_mirror_length=float(speed/(charge*b_centre/mass)/z_peak),
                analytic_loss_fraction=float(1-np.sqrt(1-1/ratio)), measured_loss_fraction=float(lost.mean()),
                predicted_loss_fraction_from_launch_pitch=float(predicted.mean()), misclassified=int(wrong.sum()), particles=count,
                misclassified_max_pitch_margin_deg=float(margin[wrong].max()) if wrong.any() else 0.,
                lost_through_ends=int(np.isin(pop.lost_surface, (3, 4)).sum()), lost_to_wall=int((pop.lost_surface == 1).sum()),
                omega_dt_max=0.1, bounce_periods=1.5)


def coil_check():
    """Quadrupole at the midplane against four infinite line currents (2D complex formula); dossier current
    for 7% of 0.35 T at the wall; Faraday's law for the RMF antenna by central differences."""
    quad = quadrupole(1000.)
    a, rw = 0.20, 0.1524
    th = np.linspace(0, 2*np.pi, 73)[:-1]
    X = np.column_stack((rw*np.cos(th), rw*np.sin(th), np.zeros_like(th)))
    b_fin = np.linalg.norm(quad.fields(X, 0.)[0], axis=1)
    zeta = rw*np.exp(1j*th)
    wires = a*np.exp(1j*(np.pi/4+np.arange(4)*np.pi/2))
    signs = np.array([1, -1, 1, -1])
    b_inf = MU0*1000/(2*np.pi)*abs((signs[None]/(zeta[:, None]-wires[None])).sum(axis=1))
    per_amp_wall = float(b_fin.mean()/1000)
    antenna, i0 = rmf_antenna(1e6, 5e-3)
    rmf = Superposition(*antenna)
    h, t, dts = 1e-4, 0.137e-6, 1e-10
    pts = np.array([[0.03, 0.02, 0.1], [-0.06, 0.04, -0.3], [0.1, -0.05, 0.5]])
    worst = 0.
    for p in pts:
        J = np.zeros((3, 3))
        for k in range(3):
            dp = np.zeros(3); dp[k] = h
            J[:, k] = (rmf.fields((p+dp)[None], t)[1][0]-rmf.fields((p-dp)[None], t)[1][0])/(2*h)
        curl = np.array([J[2, 1]-J[1, 2], J[0, 2]-J[2, 0], J[1, 0]-J[0, 1]])
        dbdt = (rmf.fields(p[None], t+dts)[0][0]-rmf.fields(p[None], t-dts)[0][0])/(2*dts)
        worst = max(worst, float(np.linalg.norm(curl+dbdt)/np.linalg.norm(dbdt)))
    return dict(quad_finite_over_infinite_mean=float(b_fin.mean()/b_inf.mean()),
                quad_azimuthal_ripple=float(b_fin.max()/b_fin.min()-1),
                quad_tesla_per_ka_at_wall=per_amp_wall*1e3,
                dossier_current_ka_for_7pct_of_035T_at_wall=float(0.07*0.35/per_amp_wall/1e3),
                rmf_current_ka_for_5mT=i0/1e3, rmf_faraday_relative_residual=worst)


def conservation_check(run, rng, count=200):
    """Energy in static fields (frozen FRC + quadrupole) and p_theta in axisymmetric fields (frozen and
    time-dependent), for omega_c,max dt = 0.2, 0.1, 0.05, 0.025."""
    t0, span = run.times[0], 4e-6
    pop0, region = run.load(count, run.ion, t0, rng)
    scale = QE*abs(region["trapped_flux_wb"])/(2*np.pi)          # e psi_O: p_theta scale
    frozen = FluxField(run.grid.r, run.grid.z, [t0], run.psi_full[0])
    b_e = run.external_field(t0)
    quad = quadrupole(QUAD_FRACTION*b_e/multipole_per_ampere(run.case.wall_radius*(1-1e-9)))
    cases = dict(static_quadrupole=(Superposition(frozen, quad), "energy"), static_axisymmetric=(frozen, "p_theta"),
                 time_dependent=(run.field, "p_theta"))
    results = {}
    for name, (fields, invariant) in cases.items():
        rows = []
        for halvings, omega_dt in enumerate((0.2, 0.1, 0.05, 0.025)):
            dt, bmax = time_step(fields, run.walls, [t0], run.ion, omega_dt, 1e-7)
            steps = int(round(span/dt))
            pop = Population(run.ion, pop0.initial["x"], pop0.initial["u"])
            if invariant == "energy":
                ref = pop.energy_ev()
                worst = [0.]

                def observe(t, p, B, ref=ref, worst=worst):
                    worst[0] = max(worst[0], float(np.max(abs(p.energy_ev()/ref[p.ids]-1))))
            else:
                ref = canonical_angular_momentum(pop, fields, t0)
                worst = [0.]

                def observe(t, p, B, ref=ref, worst=worst, fields=fields):
                    if round((t-t0)/dt) % 4 == 0:
                        worst[0] = max(worst[0], float(np.max(abs(canonical_angular_momentum(p, fields, t)-ref[p.ids]))/scale))
            track(pop, fields, t0, dt, steps, walls=run.walls, observe=observe)
            rows.append(dict(omega_dt=omega_dt, dt_ns=dt*1e9, steps=steps, survivors=int(len(pop.ids)), max_error=worst[0]))
        for a, b in zip(rows, rows[1:]):
            b["order"] = float(np.log(a["max_error"]/b["max_error"])/np.log(a["dt_ns"]/b["dt_ns"])) if invariant == "p_theta" else None
        results[name] = dict(invariant=invariant, normalisation="relative" if invariant == "energy" else "e psi_O",
                             window_us=[t0*1e6, (t0+span)*1e6], rows=rows)
    return results


def field_check(run):
    """Spline B against the display's finite-difference B (np.gradient) on the interior gas nodes, and the
    Faraday E_theta (central difference of the saved flux about saved state 1) against the ideal-MHD
    -(v x B)_theta of that state's velocity, at cell centres inside the separatrix."""
    t = run.times[0]
    br, bz = grid_fields(run.grid, run.psi_full[0])
    c = run.case
    rows = slice(run.i0+1, run.i0+c.nr)                 # interior gas rows (axis excluded)
    R, Z = np.meshgrid(run.grid.r[rows], run.grid.z[2:-2], indexing="ij")
    X = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    B = run.field.fields(X, t)[0]
    ref = np.column_stack((br[rows, 2:-2].ravel(), bz[rows, 2:-2].ravel()))
    delta = np.column_stack((B[:, 0], B[:, 2]))-ref
    scale = float(np.abs(ref).max())
    k, psi, rho, p = run.frames[1]
    t1, h = run.times[1], 1e-12
    region = run.region(1)
    rc = .5*(run.rf[:-1]+run.rf[1:]); zc = .5*(run.zf[:-1]+run.zf[1:])
    psi_c = .25*(psi[:-1, :-1]+psi[1:, :-1]+psi[:-1, 1:]+psi[1:, 1:])
    inside = psi_c < region["separatrix_psi"]
    inside[:, :region["lo"]] = False; inside[:, region["hi"]:] = False
    i, j = np.nonzero(inside)
    Xc = np.column_stack((rc[i], np.zeros(i.size), zc[j]))
    Bc = run.field.fields(Xc, t1)[0]
    e_theta = .5*(run.field.fields(Xc, t1-h)[1][:, 1]+run.field.fields(Xc, t1+h)[1][:, 1])   # phi = 0: theta_hat = y_hat
    vr, vz = run.velocity[k]
    ideal = vr[i, j]*Bc[:, 2]-vz[i, j]*Bc[:, 0]          # -(v x B)_theta = v_r B_z - v_z B_r
    rms = lambda a: float(np.sqrt(np.mean(a**2)))
    return dict(spline_vs_finite_difference_B_max_relative=float(np.abs(delta).max()/scale),
                spline_vs_finite_difference_B_rms_relative=rms(delta)/scale, gas_flux_reproduction=run.reproduction,
                faraday_rms_v_per_m=rms(e_theta), ideal_rms_v_per_m=rms(ideal),
                faraday_minus_ideal_rms_over_faraday_rms=rms(e_theta-ideal)/rms(e_theta),
                correlation=float(np.corrcoef(e_theta, ideal)[0, 1]), cells=int(i.size), time_us=t1*1e6)


def verify(out, run="outputs/fine_uniform_h2"):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    start = time.time()
    result = dict(gyration=gyration_check()); print(json.dumps(result["gyration"], indent=1), flush=True)
    result["drifts"] = drift_check(); print(json.dumps(result["drifts"], indent=1), flush=True)
    result["mirror"] = mirror_check(rng); print(json.dumps(result["mirror"], indent=1), flush=True)
    result["coils"] = coil_check(); print(json.dumps(result["coils"], indent=1), flush=True)
    frc = Run(run, 10e-6, 14.5e-6)
    result["fields"] = field_check(frc); print(json.dumps(result["fields"], indent=1), flush=True)
    result["conservation"] = conservation_check(frc, rng); print(json.dumps(result["conservation"], indent=1), flush=True)
    result["wall_s"] = time.time()-start
    (out/"verification.json").write_text(json.dumps(result, indent=1, default=float))
    return result


# ---------------------------------------------------------------- ion science runs

class OrbitRecorder:
    """Per-particle orbit diagnostics during the classification window and decimated trails afterwards."""

    def __init__(self, n, window_end, trail_ids, trail_every, steps, t0, dt):
        self.window_end, self.t0, self.dt = window_end, t0, dt
        self.bz_sign = np.zeros(n); self.null_crossings = np.zeros(n, int)
        self.lz_sign = np.zeros(n); self.lz_flips = np.zeros(n, int)
        self.mu_min = np.full(n, np.inf); self.mu_max = np.zeros(n)
        self.trail_ids = np.asarray(trail_ids)
        self.trail_slot = {int(i): s for s, i in enumerate(self.trail_ids)}
        self.every = trail_every
        self.trails = np.full((len(self.trail_ids), steps//trail_every+1, 3), np.nan, dtype=np.float32)
        self.step = 0

    def start(self, pop):
        for s, i in enumerate(self.trail_ids):
            self.trails[s, 0] = pop.x[i]

    def __call__(self, t, p, B):
        self.step += 1
        ids = p.ids
        if t <= self.window_end+1e-15:
            sb = np.sign(B[:, 2])
            prev = self.bz_sign[ids]
            self.null_crossings[ids] += (prev != 0) & (sb != prev)
            self.bz_sign[ids] = sb
            lz = np.sign(p.x[:, 0]*p.u[:, 1]-p.x[:, 1]*p.u[:, 0])
            prev = self.lz_sign[ids]
            self.lz_flips[ids] += (prev != 0) & (lz != prev)
            self.lz_sign[ids] = lz
            b = np.linalg.norm(B, axis=1)
            upar = np.einsum("nk,nk->n", p.u, B)/b
            mu = (np.einsum("nk,nk->n", p.u, p.u)-upar**2)/(2*b)
            self.mu_min[ids] = np.minimum(self.mu_min[ids], mu)
            self.mu_max[ids] = np.maximum(self.mu_max[ids], mu)
        if self.step % self.every == 0:
            sel = np.isin(ids, self.trail_ids)
            for i, x in zip(ids[sel], p.x[sel]):
                self.trails[self.trail_slot[int(i)], self.step//self.every] = x

    def classify(self, pop, mu_threshold=1.5):
        lost_early = pop.lost_time <= self.window_end
        cross = self.null_crossings >= 2
        adiabatic = self.mu_max/np.maximum(self.mu_min, 1e-300) <= mu_threshold
        cls = np.where(cross & (self.lz_flips == 0), 0, np.where(cross, 1, np.where(adiabatic, 3, 2)))
        cls[lost_early] = 4
        return cls


def run_population(run, pop0, fields, t0, t1, dt, window_s, trail_count, trail_dt):
    pop = Population(pop0.species, pop0.initial["x"], pop0.initial["u"])
    steps = int(round((t1-t0)/dt))
    every = max(1, int(round(trail_dt/dt)))
    rec = OrbitRecorder(pop.count, t0+window_s, np.arange(min(trail_count, pop.count)), every, steps, t0, dt)
    rec.start(pop)
    start = time.time()
    track(pop, fields, t0, dt, steps, walls=run.walls, observe=rec)
    return pop, rec, time.time()-start, steps


def s_bar(run, region, t, species):
    """s-bar three ways: repository convention r_s/(4 rho_i,e) (frc_scaling.py), the plain ratio
    r_s/rho_i,e, and the integral int_R^r_s r dr/(r_s rho_i(r)) with the local midplane |B| and T_i."""
    mass = SPECIES[species][0]
    ti = region["total_temperature_ev"]/2
    b_e = run.external_field(t)
    rho_e = np.sqrt(2*mass*ti*QE)/(QE*b_e)
    rs = region["separatrix_mid_m"]
    r = np.linspace(1e-4, rs, 400)
    X = np.column_stack((r, np.zeros_like(r), np.zeros_like(r)))
    b = np.linalg.norm(run.field.fields(X, t)[0], axis=1)
    rho = np.sqrt(2*mass*ti*QE)/(QE*b)
    integrand = r/(rs*rho)
    null = r[np.argmin(b)]
    outer = r >= null
    integral = float(np.sum(.5*(integrand[outer][1:]+integrand[outer][:-1])*np.diff(r[outer])))
    full = float(np.sum(.5*(integrand[1:]+integrand[:-1])*np.diff(r)))
    return dict(ti_ev=ti, b_e_t=b_e, rho_i_external_m=rho_e, separatrix_mid_m=rs, null_radius_m=float(null),
                s_bar_repository=rs/(4*rho_e), rs_over_rho_i=rs/rho_e, s_bar_local_null_to_separatrix=integral,
                s_bar_local_axis_to_separatrix=full)


def confinement_curve(pop, t0, t1, classes=None, points=201):
    t = np.linspace(t0, t1, points)
    lost = np.where(np.isfinite(pop.lost_time), pop.lost_time, np.inf)
    out = dict(time_us=(t*1e6).round(3).tolist(), all=[float((lost > s).mean()) for s in t])
    if classes is not None:
        for c, name in enumerate(CLASSES):
            m = classes == c
            if m.any():
                out[name] = [float((lost[m] > s).mean()) for s in t]
    return out


def loss_summary(pop):
    lost = np.isfinite(pop.lost_time)
    out = dict(lost_fraction=float(lost.mean()))
    for code in (1, 2, 3, 4):
        m = pop.lost_surface == code
        if m.any():
            out[SURFACES[code]] = dict(fraction=float(m.mean()), mean_energy_ev=float(pop.lost_energy_ev[m].mean()),
                                       median_time_us=float(np.median(pop.lost_time[m])*1e6))
    return out


def sector_amplitude(pop, m=2):
    """|<exp(i m phi)>| of the loss azimuths, separately for wall losses at z > 0 and z < 0 and for each end:
    the m = 2 imprint of a quadrupole, m = 4 of an octupole (the field-line divergence reverses sign between the
    two halves, i.e. the pattern is rotated by 180/m deg). noise = 1/sqrt(count) for uniformly spread azimuths."""
    out = {}
    z = pop.lost_position[:, 2]
    for label, sel in (("wall z>0", (pop.lost_surface == 1) & (z > 0)), ("wall z<0", (pop.lost_surface == 1) & (z < 0)),
                       ("end +z", pop.lost_surface == 3), ("end -z", pop.lost_surface == 4)):
        if sel.sum() >= 10:
            phi = np.arctan2(pop.lost_position[sel, 1], pop.lost_position[sel, 0])
            mean = np.exp(1j*m*phi).mean()
            out[label] = dict(count=int(sel.sum()), amplitude=float(abs(mean)), phase_deg=float(np.degrees(np.angle(mean))/m),
                              noise=float(1/np.sqrt(sel.sum())))
    return out


def ions(run_path, out, start_us=10., end_us=30., count=4000, seed=1, window_us=2., trails=300, omega_dt=0.1):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    t0, t1 = start_us*1e-6, end_us*1e-6
    run = Run(run_path, t0, t1)
    t1 = min(t1, run.times[-1])
    rng = np.random.default_rng(seed)
    pop0, region = run.load(count, run.ion, t0, rng)
    b_e = run.external_field(t0)
    quad_currents = {"wall": QUAD_FRACTION*b_e/multipole_per_ampere(run.case.wall_radius*(1-1e-9)),
                     "separatrix": QUAD_FRACTION*b_e/multipole_per_ampere(region["separatrix_mid_m"])}
    cases = {"none": run.field}
    for where, current in quad_currents.items():
        cases[f"quadrupole 7% B_e at the {where}"] = Superposition(run.field, quadrupole(current))
    sbar = s_bar(run, region, t0, run.ion)
    summary = dict(run=str(run_path), species=run.ion, count=count, start_us=start_us, end_us=t1*1e6, seed=seed,
                   classification_window_us=window_us, region={k: region[k] for k in ("separatrix_mid_m", "length_m",
                   "trapped_flux_wb", "density_m3", "total_temperature_ev", "separatrix_psi")},
                   b_e_t=b_e, quadrupole_bar_current_a=quad_currents, s_bar=sbar,
                   flux_reproduction=run.reproduction, cases={})
    print(json.dumps({k: summary[k] for k in ("run", "b_e_t", "quadrupole_bar_current_a", "s_bar")}, default=float), flush=True)
    results = {}
    for name, fields in cases.items():
        dt, bmax = time_step(fields, run.walls, [t0, t1], run.ion, omega_dt, 1e-7)
        pop, rec, wall_s, steps = run_population(run, pop0, fields, t0, t1, dt, window_us*1e-6, trails, 20e-9)
        classes = rec.classify(pop)
        results[name] = (pop, rec, classes)
        np.savez_compressed(out/f"ions_{'none' if name == 'none' else name.split()[-1]}.npz", x0=pop.initial["x"], u0=pop.initial["u"],
                            lost_time=pop.lost_time, lost_position=pop.lost_position, lost_surface=pop.lost_surface,
                            lost_energy_ev=pop.lost_energy_ev, orbit_class=classes, null_crossings=rec.null_crossings,
                            lz_flips=rec.lz_flips, mu_min=rec.mu_min, mu_max=rec.mu_max)
        frac = {n: float((classes == c).mean()) for c, n in enumerate(CLASSES)}
        alive = ~np.isfinite(pop.lost_time)
        summary["cases"][name] = dict(dt_ns=dt*1e9, steps=steps, b_max_t=bmax, wall_s=wall_s, class_fractions=frac,
                                      confined_at_end=float(alive.mean()),
                                      confined_by_class={n: float(alive[classes == c].mean()) for c, n in enumerate(CLASSES) if (classes == c).any()},
                                      losses=loss_summary(pop), loss_m2=sector_amplitude(pop),
                                      confinement=confinement_curve(pop, t0, t1, classes),
                                      mu_ratio_quantiles=np.quantile(rec.mu_max/np.maximum(rec.mu_min, 1e-300), [.1, .25, .5, .75, .9]).tolist(),
                                      class_fractions_mu_threshold_1p2=_fractions(rec.classify(pop, 1.2)),
                                      class_fractions_mu_threshold_2=_fractions(rec.classify(pop, 2.)))
        print(name, json.dumps({k: v for k, v in summary["cases"][name].items() if k != "confinement"}, default=float), flush=True)
    base = results["none"][0]
    for name in list(cases)[1:]:
        pop = results[name][0]
        a, b = ~np.isfinite(base.lost_time), ~np.isfinite(pop.lost_time)
        summary["cases"][name]["paired"] = dict(confined_without_lost_with=float((a & ~b).mean()),
                                                lost_without_confined_with=float((~a & b).mean()),
                                                change_in_confined_fraction=float(b.mean()-a.mean()))
    (out/"summary.json").write_text(json.dumps(summary, indent=1, default=float))
    figures(out, run, summary, results, t0, t1)
    write_trails(run, results, t0, out)
    return summary


def _fractions(classes):
    return {n: float((classes == c).mean()) for c, n in enumerate(CLASSES)}


def figures(out, run, summary, results, t0, t1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(results)
    fig, axes = plt.subplots(1, len(names), figsize=(5*len(names), 3.6), sharey=True)
    for ax, name in zip(np.atleast_1d(axes), names):
        curve = summary["cases"][name]["confinement"]
        ax.plot(curve["time_us"], curve["all"], color="k", lw=2, label="all")
        for cls in CLASSES:
            if cls in curve:
                ax.plot(curve["time_us"], curve[cls], lw=1.2, label=cls, color={"betatron": "#3987e5", "mixed": "#8a5cc2", "drift": "#d55181",
                                                                            "cyclotron": "#199e70", "lost in window": "#c98500"}[cls])
        ax.set_title(name, fontsize=9); ax.set_xlabel("time (µs)"); ax.grid(alpha=.3)
    np.atleast_1d(axes)[0].set_ylabel("confined fraction"); np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out/"confinement.png", dpi=130); plt.close(fig)
    fig, axes = plt.subplots(2, len(names), figsize=(5*len(names), 7))
    axes = np.array(axes).reshape(2, -1)
    for col, name in enumerate(names):
        pop = results[name][0]
        wall = pop.lost_surface == 1
        ends = np.isin(pop.lost_surface, (3, 4))
        ax = axes[0, col]
        ax.scatter(pop.lost_position[wall, 2], np.degrees(np.arctan2(pop.lost_position[wall, 1], pop.lost_position[wall, 0])), s=3, c=pop.lost_time[wall]*1e6, cmap="viridis")
        ax.set_xlim(-run.walls.z_end, run.walls.z_end); ax.set_ylim(-180, 180); ax.set_xlabel("z (m)"); ax.set_ylabel("azimuth (deg)")
        ax.set_title(f"{name}\nwall losses: {wall.sum()}", fontsize=9)
        ax = axes[1, col]
        for code, marker in ((3, "o"), (4, "x")):
            m = pop.lost_surface == code
            ax.scatter(pop.lost_position[m, 0]*100, pop.lost_position[m, 1]*100, s=4, marker=marker, c=pop.lost_time[m]*1e6, cmap="viridis",
                       label=SURFACES[code])
        circle = np.linspace(0, 2*np.pi, 200)
        ax.plot(run.walls.r_outer*100*np.cos(circle), run.walls.r_outer*100*np.sin(circle), color="#888", lw=.8)
        if run.walls.r_inner > 0:
            ax.plot(run.walls.r_inner*100*np.cos(circle), run.walls.r_inner*100*np.sin(circle), color="#888", lw=.8)
        ax.set_aspect("equal"); ax.set_xlabel("x (cm)"); ax.set_ylabel("y (cm)")
        ax.set_title(f"end losses: {ends.sum()} (o +z, x −z), colour = time", fontsize=9)
    fig.tight_layout(); fig.savefig(out/"loss_maps.png", dpi=130); plt.close(fig)
    pop, rec, classes = results["none"]
    x0 = pop.initial["x"]
    fig, ax = plt.subplots(figsize=(8, 3.4))
    for c, cls in enumerate(CLASSES):
        m = classes == c
        ax.scatter(x0[m, 2], np.hypot(x0[m, 0], x0[m, 1])*100, s=2, label=f"{cls} {m.mean():.0%}")
    ax.set_xlabel("z (m)"); ax.set_ylabel("r (cm)"); ax.legend(fontsize=7, markerscale=4); ax.set_title("Initial positions by orbit class (no quadrupole)")
    fig.tight_layout(); fig.savefig(out/"orbit_classes.png", dpi=130); plt.close(fig)


def write_trails(run, results, t0, out, per_class=5):
    """RUN/particles.json: decimated actual orbits (positions every 20 ns, metres) for the viewer, a few per
    orbit class and a few of the lost ions, for each case."""
    sets = []
    for name, (pop, rec, classes) in results.items():
        chosen = []
        for c, cls in enumerate(CLASSES[:4]):
            ids = [i for i in rec.trail_ids if classes[i] == c and not np.isfinite(pop.lost_time[i])][:per_class]
            chosen += ids
        lost = [i for i in rec.trail_ids if np.isfinite(pop.lost_time[i]) and i not in chosen][:per_class]
        trails = []
        for i in chosen+lost:
            pts = rec.trails[rec.trail_slot[int(i)]]
            n = int(np.isfinite(pts[:, 0]).sum())
            pts = pts[:n].astype(float)
            if np.isfinite(pop.lost_time[i]):
                pts = np.vstack((pts, pop.lost_position[i]))
            trails.append(dict(id=int(i), **{"class": CLASSES[classes[i]]},
                               lost_us=None if not np.isfinite(pop.lost_time[i]) else round(float(pop.lost_time[i]*1e6), 3),
                               surface=SURFACES[pop.lost_surface[i]], points=np.round(pts, 4).tolist()))
        sets.append(dict(label="No quadrupole" if name == "none" else name.replace("quadrupole", "Quadrupole"), trails=trails))
    data = dict(source="examples/particle_tracking.py", run=run.path.name, species=run.ion, t0_us=round(t0*1e6, 6),
                sample_us=rec.every*rec.dt*1e6,
                colors=COLORS, sets=sets,
                note="Actual integrated full orbits (relativistic Boris) of thermal ions in the saved fields; positions every "
                     f"{rec.every*rec.dt*1e9:.0f} ns. Test particles: no self-consistent field response.")
    text = json.dumps(data, separators=(",", ":"))
    (run.path/"particles.json").write_text(text)
    (out/"particles.json").write_text(text)


# ---------------------------------------------------------------- multipole strength/order scan

ORDERS = {"q": 2, "o": 4}        # case codes: q7 = 4-bar quadrupole at 7% B_e, o7 = 8-bar octupole at 7% B_e


def _multipole_case(field, walls, species, x0, u0, order, current, t0, t1, omega_dt):
    """One scan case (worker process): the FRC field alone (order 0) or plus a static 2*order-bar multipole."""
    fields = field if order == 0 else Superposition(field, multipole(current, order))
    dt, bmax = time_step(fields, walls, [t0, t1], species, omega_dt, 1e-7)
    pop = Population(species, x0, u0)
    steps = int(round((t1-t0)/dt))
    start = time.time()
    track(pop, fields, t0, dt, steps, walls=walls)
    return dict(dt=dt, b_max=bmax, steps=steps, wall_s=time.time()-start, lost_time=pop.lost_time,
                lost_position=pop.lost_position, lost_surface=pop.lost_surface, lost_energy_ev=pop.lost_energy_ev)


def multipoles(run_path, out, start_us=10., end_us=30., count=4000, seed=1, cases="none,q3,q5,q7,q9,o7", omega_dt=0.1,
               workers=0):
    """Ion-loss cost of n=2 stabilization against multipole strength and order. The same ions (seed) are tracked
    in the FRC field alone and with each static multipole, whose bar current sets the mean |B| around the wall
    circle at the midplane to the given percentage of B_e. Vacuum superposition: no screening by the plasma."""
    from concurrent.futures import ProcessPoolExecutor
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    t0, t1 = start_us*1e-6, end_us*1e-6
    run = Run(run_path, t0, t1)
    t1 = min(t1, run.times[-1])
    pop0, region = run.load(count, run.ion, t0, np.random.default_rng(seed))
    b_e = run.external_field(t0)
    wall, rs = run.case.wall_radius*(1-1e-9), region["separatrix_mid_m"]
    spec = {}
    for name in [c.strip() for c in str(cases).split(",") if c.strip()]:
        if name == "none":
            spec[name] = dict(order=0, percent_b_e=0., bar_current_a=0.)
            continue
        order, percent = ORDERS[name[0]], float(name[1:])
        current = percent/100*b_e/multipole_per_ampere(wall, order)
        spec[name] = dict(order=order, bars=2*order, percent_b_e=percent, bar_current_a=current,
                          wall_max_over_b_e=current*multipole_per_ampere(wall, order, np.max)/b_e,
                          separatrix_mean_over_b_e=current*multipole_per_ampere(rs, order)/b_e)
    print(json.dumps(dict(run=str(run_path), b_e_t=b_e, separatrix_mid_m=rs, cases=spec), default=float), flush=True)
    with ProcessPoolExecutor(max_workers=workers or len(spec)) as pool:
        jobs = {name: pool.submit(_multipole_case, run.field, run.walls, run.ion, pop0.initial["x"], pop0.initial["u"],
                                  c["order"], c["bar_current_a"], t0, t1, omega_dt) for name, c in spec.items()}
        results = {name: job.result() for name, job in jobs.items()}
    summary = dict(run=str(run_path), species=run.ion, count=count, start_us=start_us, end_us=t1*1e6, seed=seed,
                   b_e_t=b_e, wall_radius_m=run.case.wall_radius, separatrix_mid_m=rs, bar_radius_m=0.20,
                   model="vacuum superposition of static Biot-Savart multipoles on the saved axisymmetric fields; "
                         "no plasma screening, no collisions, test particles", cases={})
    pops = {}
    for name, res in results.items():
        pop = pops[name] = Population(run.ion, pop0.initial["x"], pop0.initial["u"])
        for key in ("lost_time", "lost_position", "lost_surface", "lost_energy_ev"):
            setattr(pop, key, res[key])
        np.savez_compressed(out/f"ions_{name}.npz", x0=pop.initial["x"], u0=pop.initial["u"], lost_time=pop.lost_time,
                            lost_position=pop.lost_position, lost_surface=pop.lost_surface, lost_energy_ev=pop.lost_energy_ev)
        lost = np.where(np.isfinite(pop.lost_time), pop.lost_time, np.inf)
        wall_hit = pop.lost_surface == 1
        summary["cases"][name] = dict(
            spec[name], dt_ns=res["dt"]*1e9, steps=res["steps"], b_max_t=res["b_max"], wall_s=res["wall_s"],
            confined={f"{t:g}": float((lost > t*1e-6).mean()) for t in [*np.arange(start_us, t1*1e6-1e-9, 5.), t1*1e6]},
            lost_to_wall=float(wall_hit.mean()), lost_to_ends=float(np.isin(pop.lost_surface, (3, 4)).mean()),
            lost_to_inner_tube=float((pop.lost_surface == 2).mean()), losses=loss_summary(pop),
            wall_loss_abs_z_quantiles_m=(np.quantile(abs(pop.lost_position[wall_hit, 2]), [.1, .5, .9]).tolist()
                                         if wall_hit.any() else None),
            sector_m2=sector_amplitude(pop, 2), sector_m4=sector_amplitude(pop, 4),
            confinement=confinement_curve(pop, t0, t1))
    if "none" in pops:
        a = ~np.isfinite(pops["none"].lost_time)
        for name, pop in pops.items():
            b = ~np.isfinite(pop.lost_time)
            gone, gained = (a & ~b).mean(), (~a & b).mean()
            summary["cases"][name]["paired"] = dict(
                confined_without_lost_with=float(gone), lost_without_confined_with=float(gained),
                change_in_confined_fraction=float(b.mean()-a.mean()),
                change_standard_error=float(np.sqrt(max(gone+gained-(gained-gone)**2, 0.)/count)))
    (out/"summary.json").write_text(json.dumps(summary, indent=1, default=float))
    print(f"{'case':6} {'I_bar kA':>8} {'sep/B_e':>7} {'conf20':>7} {'conf_end':>8} {'d_end':>13} {'wall':>6} {'ends':>6} "
          f"{'m2 z>0/z<0':>11} {'m4 z>0/z<0':>11}")
    for name, c in summary["cases"].items():
        amp = lambda key: "/".join(f"{c[key][h]['amplitude']:.2f}" if h in c[key] else " -- " for h in ("wall z>0", "wall z<0"))
        d = c.get("paired", {})
        print(f"{name:6} {c['bar_current_a']/1e3:8.2f} {c.get('separatrix_mean_over_b_e', 0):7.3f} {c['confined'].get('20', np.nan):7.3f} "
              f"{c['confined'][f'{t1*1e6:g}']:8.3f} {d.get('change_in_confined_fraction', 0):+7.3f}±{d.get('change_standard_error', 0):.3f} "
              f"{c['lost_to_wall']:6.3f} {c['lost_to_ends']:6.3f} {amp('sector_m2'):>11} {amp('sector_m4'):>11}", flush=True)
    multipole_figure(out, summary, pops, run)
    return summary


def multipole_figure(out, summary, pops, run):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    colors = {"none": "#222222", "q3": "#9ec5f0", "q5": "#5b9be0", "q7": "#2a6cc0", "q9": "#123f80", "o7": "#d55181"}
    for name, c in summary["cases"].items():
        curve = c["confinement"]
        axes[0].plot(curve["time_us"], curve["all"], color=colors.get(name), lw=1.6, label=name)
    axes[0].set_xlabel("time (µs)"); axes[0].set_ylabel("confined fraction"); axes[0].grid(alpha=.3); axes[0].legend(fontsize=8)
    axes[0].set_title("Same ions with each multipole (q = 4-bar, o = 8-bar; % of B_e at the wall)", fontsize=9)
    end = f"{summary['end_us']:g}"
    for order, marker in ((2, "o"), (4, "s")):
        pts = sorted((c["percent_b_e"], c["confined"][end]) for c in summary["cases"].values() if c["order"] in (0, order))
        axes[1].plot(*zip(*pts), marker=marker, label={2: "quadrupole", 4: "octupole"}[order])
    axes[1].set_xlabel("mean |B| at the wall (% of B_e)"); axes[1].set_ylabel(f"confined at {end} µs"); axes[1].grid(alpha=.3)
    axes[1].legend(fontsize=8)
    bins = np.linspace(-180, 180, 73)
    for name, pop in pops.items():
        if name == "none":
            continue
        for half, style in ((1, "-"), (-1, "--")):
            sel = (pop.lost_surface == 1) & (half*pop.lost_position[:, 2] > 0)
            phi = np.degrees(np.arctan2(pop.lost_position[sel, 1], pop.lost_position[sel, 0]))
            h, _ = np.histogram(phi, bins)
            axes[2].step(bins[:-1], h, where="post", ls=style, color=colors.get(name), lw=1, label=f"{name} {'z>0' if half > 0 else 'z<0'}")
    axes[2].set_xlabel("azimuth of wall loss (deg)"); axes[2].set_ylabel("ions per 5°"); axes[2].legend(fontsize=6, ncol=2)
    axes[2].set_title("Wall-loss sectors (vacuum superposition: upper bound)", fontsize=9)
    fig.tight_layout(); fig.savefig(out/"multipoles.png", dpi=130); plt.close(fig)


# ---------------------------------------------------------------- RMF electrons

def rmf(run_path, out, start_us=10., duration_us=2., count=300, frequency=1e6, amplitude=5e-3, parity="even",
        radius=0.16, half_length=0.6, omega_dt=0.3, seed=3):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    t0, t1 = start_us*1e-6, (start_us+duration_us)*1e-6
    run = Run(run_path, t0, t1+0.1e-6)
    rng = np.random.default_rng(seed)
    pop0, region = run.load(count, "e-", t0, rng)
    antenna, i0 = rmf_antenna(frequency, amplitude, parity, radius, half_length)
    omega = 2*np.pi*frequency
    # where is omega above the local ion cyclotron frequency? (fraction of loaded particle positions)
    b_local = np.linalg.norm(run.field.fields(pop0.initial["x"], t0)[0], axis=1)
    wci = QE*b_local/SPECIES[run.ion][0]
    summary = dict(run=str(run_path), frequency_hz=frequency, amplitude_t=amplitude, parity=parity, antenna_radius_m=radius,
                   antenna_half_length_m=half_length, antenna_current_a=i0, count=count, start_us=start_us, duration_us=duration_us,
                   fraction_of_loaded_positions_with_omega_above_omega_ci=float((omega > wci).mean()),
                   ion_cyclotron_frequency_at_b_e_hz=float(QE*run.external_field(t0)/SPECIES[run.ion][0]/(2*np.pi)),
                   electron_cyclotron_frequency_at_b_e_hz=float(QE*run.external_field(t0)/ME/(2*np.pi)), cases={})
    per = {}
    for name, fields in (("no RMF", run.field), ("RMF", Superposition(run.field, *antenna))):
        dt, bmax = time_step(fields, run.walls, [t0], "e-", omega_dt, 1e-7)
        steps = int(round((t1-t0)/dt))
        pop = Population("e-", pop0.initial["x"], pop0.initial["u"])
        sample = max(1, int(round(2e-9/dt)))
        rec = dict(theta_dot=np.full((count, steps//sample+1), np.nan, np.float32), r=np.full((count, steps//sample+1), np.nan, np.float32),
                   bz_over_bw=np.full((count, steps//sample+1), np.nan, np.float32), k=0, n=0)

        def observe(t, p, B, rec=rec):
            rec["k"] += 1
            if rec["k"] % sample:
                return
            rec["n"] += 1
            x, u = p.x, p.u/gamma_of(p.u)[:, None]
            r2 = x[:, 0]**2+x[:, 1]**2
            rec["theta_dot"][p.ids, rec["n"]] = (x[:, 0]*u[:, 1]-x[:, 1]*u[:, 0])/r2
            rec["r"][p.ids, rec["n"]] = np.sqrt(r2)
            rec["bz_over_bw"][p.ids, rec["n"]] = np.abs(B[:, 2])/amplitude
        start = time.time()
        track(pop, fields, t0, dt, steps, walls=run.walls, observe=observe)
        settle = int(rec["n"]*0.25)                         # skip the first quarter of the run
        td, rr, bz = (rec[k][:, settle:rec["n"]+1] for k in ("theta_dot", "r", "bz_over_bw"))
        mean_rate = np.nanmean(td, axis=1)/omega           # per-electron time-averaged theta-dot / omega
        weak = np.nanmean(bz, axis=1) < 1.                 # time-averaged |B_z| below the RMF amplitude
        bins = np.linspace(0, run.case.wall_radius, 9)
        rmean = np.nanmean(rr, axis=1)
        profile = []
        for a, b in zip(bins[:-1], bins[1:]):
            m = (rmean >= a) & (rmean < b) & np.isfinite(mean_rate)
            if m.sum() >= 3:
                profile.append(dict(r_cm=[round(a*100, 2), round(b*100, 2)], electrons=int(m.sum()),
                                    mean_theta_dot_over_omega=float(np.mean(mean_rate[m])),
                                    median_theta_dot_over_omega=float(np.median(mean_rate[m]))))
        alive = ~np.isfinite(pop.lost_time)
        energy = np.full(count, np.nan); energy[pop.ids] = pop.energy_ev()
        per[name] = dict(rate=mean_rate, r=rmean, bz_over_bw=np.nanmean(bz, axis=1), energy_ev=energy)
        np.savez_compressed(out/f"rmf_{name.replace(' ', '_')}.npz", **per[name])
        summary["cases"][name] = dict(dt_ps=dt*1e12, steps=steps, b_max_t=bmax, wall_s=time.time()-start,
                                      lost_fraction=float(1-alive.mean()), losses=loss_summary(pop),
                                      mean_energy_ev_start=float(Population("e-", pop0.initial["x"], pop0.initial["u"]).energy_ev().mean()),
                                      mean_energy_ev_end=float(pop.energy_ev().mean()) if len(pop.ids) else None,
                                      median_energy_ev_end=float(np.median(pop.energy_ev())) if len(pop.ids) else None,
                                      theta_dot_over_omega_all=dict(mean=float(np.nanmean(mean_rate)), median=float(np.nanmedian(mean_rate))),
                                      theta_dot_over_omega_weak_bz=dict(count=int((weak & np.isfinite(mean_rate)).sum()),
                                                                        median=float(np.nanmedian(mean_rate[weak])) if weak.any() else None),
                                      theta_dot_over_omega_strong_bz=dict(count=int((~weak & np.isfinite(mean_rate)).sum()),
                                                                          median=float(np.nanmedian(mean_rate[~weak])) if (~weak).any() else None),
                                      entrained_fraction_within_30pct=float(np.mean(np.abs(mean_rate[np.isfinite(mean_rate)]-1) < 0.3)),
                                      radial_profile=profile)
        print(name, json.dumps(summary["cases"][name], default=float), flush=True)
    a, b = per["no RMF"], per["RMF"]
    shift = b["rate"]-a["rate"]                            # paired: same electrons with and without the RMF
    ok = np.isfinite(shift)
    weak = a["bz_over_bw"] < 1.
    summary["paired_theta_dot_shift_over_omega"] = dict(
        median_all=float(np.median(shift[ok])), median_weak_bz=float(np.median(shift[ok & weak])) if (ok & weak).any() else None,
        median_strong_bz=float(np.median(shift[ok & ~weak])), quartiles_all=np.quantile(shift[ok], [.25, .75]).tolist(),
        weak_bz_count=int((ok & weak).sum()))
    print(json.dumps(summary["paired_theta_dot_shift_over_omega"]), flush=True)
    (out/"rmf_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


def _value(text):
    for kind in (int, float):
        try:
            return kind(text)
        except ValueError:
            pass
    return text


if __name__ == "__main__":
    command, *args = sys.argv[1:]
    positional = [a for a in args if "=" not in a]
    options = {k: _value(v) for k, v in (a.split("=", 1) for a in args if "=" in a)}
    if command == "verify":
        verify(*positional, **options)
    elif command == "ions":
        ions(*positional, **options)
    elif command == "rmf":
        rmf(*positional, **options)
    elif command == "multipoles":
        multipoles(*positional, **options)
    else:
        raise SystemExit(__doc__)
