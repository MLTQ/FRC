"""Full-orbit test particles in saved axisymmetric FRC fields plus analytic 3D coil fields.

Fields (companion: particle_tracking.md):
- `FluxField`: B = grad(psi) x theta/r and E_theta = -(dpsi/dt)/r from a poloidal flux psi(r, z) saved at
  times t_k. psi = r^2 g(r, z, t) with g a cubic B-spline interpolant on the uniform (r, z) grid (even in r),
  linear in t between saved states. B is then exactly divergence-free, exactly the curl of A_theta = r g,
  and E_theta is exactly Faraday's law for that psi, so p_theta = gamma m r v_theta + q psi is an exact
  invariant of the interpolated axisymmetric fields.
- `Filaments`: closed polylines of straight current filaments with I(t). B by the finite-segment
  Biot-Savart formula, E = -dI/dt A with the closed-form segment vector potential (inductive vacuum field,
  no retardation, no plasma response).
- `Superposition`: sum of field parts.

Integrator: relativistic Boris in drift-kick-drift form (positions and momenta at the same times). Momentum
u = gamma v; a pure magnetic kick rotates u exactly, so energy is conserved to round-off in static B.

Losses: `Walls` (outer wall, optional inner tube, tube ends); the crossing point and time on the surface are
interpolated within the step.
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.constants import c as C0, e as QE, m_e as ME, m_p as MP, mu_0 as MU0, physical_constants
from scipy.ndimage import spline_filter

MD = physical_constants["deuteron mass"][0]
SPECIES = {"H+": (MP, QE), "D+": (MD, QE), "e-": (ME, -QE)}
SURFACES = ("", "wall", "inner tube", "end +z", "end -z")


def bspline_weights(s):
    """Cubic B-spline weights and their s-derivatives for the nodes floor(s)-1 .. floor(s)+2."""
    i = np.floor(s).astype(np.int64)
    t = s-i
    t2, t3, u = t*t, t*t*t, 1-t
    w = np.stack((u*u*u/6, (3*t3-6*t2+4)/6, (-3*t3+3*t2+3*t+1)/6, t3/6), axis=-1)
    dw = np.stack((-u*u/2, (3*t2-4*t)/2, (-3*t2+2*t+1)/2, t2/2), axis=-1)
    return i, w, dw


class FluxField:
    """Axisymmetric fields from psi[k] (shape (nt, nr, nz)) saved at times[k] on nodes r = 0, dr, .. and z."""

    def __init__(self, r, z, times, psi):
        r, z, psi = np.asarray(r, float), np.asarray(z, float), np.asarray(psi, float)
        if psi.ndim == 2:
            psi = psi[None]
        self.dr, self.dz, self.z0 = r[1]-r[0], z[1]-z[0], z[0]
        if abs(r[0]) > 1e-9*self.dr or not np.allclose(np.diff(r), self.dr) or not np.allclose(np.diff(z), self.dz):
            raise ValueError("FluxField needs a uniform grid with r[0] = 0")
        g = np.empty_like(psi)
        g[:, 1:] = psi[:, 1:]/r[1:, None]**2
        g[:, 0] = (16*psi[:, 1]-psi[:, 2])/(12*self.dr**2)        # psi = a r^2 + b r^4 near the axis
        coef = np.stack([spline_filter(gk, order=3, mode="mirror") for gk in g])
        self.coef = np.pad(coef, ((0, 0), (2, 2), (2, 2)), mode="reflect")   # padded index i+2 holds node i
        self.times = np.atleast_1d(np.asarray(times, float))
        self.shape = r.size, z.size
        self.r_max, self.z_lo, self.z_hi = r[-1], z[0], z[-1]
        self._cache = (None, None)

    def _blend(self, t):
        """Spline coefficients of g + i dg/dt at time t (linear in t between saved states)."""
        if self._cache[0] == t:
            return self._cache[1]
        if len(self.times) == 1:
            c = self.coef[0].astype(complex)
        else:
            if t < self.times[0]-1e-15 or t > self.times[-1]+1e-15:
                raise ValueError(f"t = {t:.4g} s is outside the saved states [{self.times[0]:.4g}, {self.times[-1]:.4g}]")
            k = int(np.clip(np.searchsorted(self.times, t, side="right")-1, 0, len(self.times)-2))
            span = self.times[k+1]-self.times[k]
            tau = (t-self.times[k])/span
            dc = self.coef[k+1]-self.coef[k]
            c = self.coef[k]+tau*dc+1j*dc/span
        self._cache = (t, c)
        return c

    def evaluate(self, x, y, z, t, derivatives=True):
        """r, g + i dg/dt and (if derivatives) dg/dr, dg/dz at the points."""
        r = np.hypot(x, y)
        nr, nz = self.shape
        i, wr, dwr = bspline_weights(np.clip(r/self.dr, 0, nr-1-1e-9))
        j, wz, dwz = bspline_weights(np.clip((z-self.z0)/self.dz, 0, nz-1-1e-9))
        c = self._blend(t)
        k = np.arange(1, 5)
        G = c[(i[:, None]+k)[:, :, None], (j[:, None]+k)[:, None, :]]          # (N, 4, 4)
        Gz = np.einsum("nj,nij->ni", wz, G)
        g = np.einsum("ni,ni->n", wr, Gz)
        if not derivatives:
            return r, g
        g_r = np.einsum("ni,ni->n", dwr, Gz.real)/self.dr
        g_z = np.einsum("ni,nj,nij->n", wr, dwz, G.real)/self.dz
        return r, g, g_r, g_z

    def fields(self, X, t):
        x, y, z = X[:, 0], X[:, 1], X[:, 2]
        r, g, g_r, g_z = self.evaluate(x, y, z, t)
        B = np.column_stack((-x*g_z, -y*g_z, 2*g.real+r*g_r))
        E = np.column_stack((y*g.imag, -x*g.imag, np.zeros_like(x)))    # E_theta = -r dg/dt
        return B, E

    def psi(self, X, t):
        r, g = self.evaluate(X[:, 0], X[:, 1], X[:, 2], t, derivatives=False)
        return r*r*g.real


class Filaments:
    """Closed polylines of straight filaments carrying I(t); `current(t)` returns (I, dI/dt)."""

    def __init__(self, paths, current, label=""):
        a = np.concatenate([np.asarray(p, float)[:-1] for p in paths])
        b = np.concatenate([np.asarray(p, float)[1:] for p in paths])
        keep = np.linalg.norm(b-a, axis=1) > 0
        self.a, self.b, self.current, self.label = a[keep], b[keep], current, label
        self.length = np.linalg.norm(self.b-self.a, axis=1)
        self.unit = (self.b-self.a)/self.length[:, None]

    def per_ampere(self, X, potential=True):
        """B and A per ampere at X (N, 3): finite-segment Biot-Savart and segment vector potential."""
        (x1, y1, z1), (x2, y2, z2) = ([X[:, k, None]-end[None, :, k] for k in range(3)] for end in (self.a, self.b))
        n1, n2 = np.sqrt(x1*x1+y1*y1+z1*z1), np.sqrt(x2*x2+y2*y2+z2*z2)
        s = n1+n2
        factor = MU0/(4*np.pi)*s/(n1*n2*(n1*n2+x1*x2+y1*y2+z1*z2))
        B = np.column_stack((np.einsum("nm,nm->n", factor, y1*z2-z1*y2), np.einsum("nm,nm->n", factor, z1*x2-x1*z2),
                             np.einsum("nm,nm->n", factor, x1*y2-y1*x2)))
        if not potential:
            return B, None
        A = MU0/(4*np.pi)*np.log((s+self.length)/(s-self.length))@self.unit
        return B, A

    def fields(self, X, t):
        current, rate = self.current(t)
        B, A = self.per_ampere(X, potential=rate != 0)
        return current*B, (-rate*A if rate != 0 else np.zeros_like(B))


class Superposition:
    """Sum of field parts; psi() from the first part that has one (the axisymmetric FRC field)."""

    def __init__(self, *parts):
        self.parts = [p for p in parts if p is not None]

    def fields(self, X, t):
        B, E = self.parts[0].fields(X, t)
        for p in self.parts[1:]:
            b, e = p.fields(X, t)
            B, E = B+b, E+e
        return B, E

    def psi(self, X, t):
        return next(p for p in self.parts if hasattr(p, "psi")).psi(X, t)


class UniformField:
    """Uniform static B and E (verification)."""

    def __init__(self, B, E=(0., 0., 0.)):
        self.B, self.E = np.asarray(B, float), np.asarray(E, float)

    def fields(self, X, t):
        return np.broadcast_to(self.B, X.shape).copy(), np.broadcast_to(self.E, X.shape).copy()


class GradientField:
    """B = B0 (1 + x/L) z_hat (verification of the grad-B drift)."""

    def __init__(self, b0, length):
        self.b0, self.length = b0, length

    def fields(self, X, t):
        B = np.zeros_like(X)
        B[:, 2] = self.b0*(1+X[:, 0]/self.length)
        return B, np.zeros_like(X)


# ---------------------------------------------------------------- coil geometry

def arc(radius, phi0, phi1, z, segments):
    a = np.linspace(phi0, phi1, segments+1)
    return np.column_stack((radius*np.cos(a), radius*np.sin(a), np.full(a.size, z)))


def rotate_z(paths, angle):
    c, s = np.cos(angle), np.sin(angle)
    m = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    return [np.asarray(p)@m.T for p in paths]


def quadrupole_paths(radius=0.20, half_length=1.1, phase=np.pi/4, arc_segments=6):
    """Four axial bars at `radius` (azimuths phase + k pi/2) in series with alternating current directions,
    joined by arcs at the bar radius at z = +/-half_length: one closed loop (BUILD_DOSSIER.md section 6)."""
    pts = []
    for k in range(4):
        phi = phase+k*np.pi/2
        z0, z1 = (-half_length, half_length) if k % 2 == 0 else (half_length, -half_length)
        bar = np.array([[radius*np.cos(phi), radius*np.sin(phi), z0], [radius*np.cos(phi), radius*np.sin(phi), z1]])
        pts.append(bar if k == 0 else bar[1:])
        pts.append(arc(radius, phi, phi+np.pi/2, z1, arc_segments)[1:])
    return [np.vstack(pts)]


def quadrupole(current, **geometry):
    """Static quadrupole: Filaments with a constant bar current [A]."""
    return Filaments(quadrupole_paths(**geometry), lambda t: (current, 0.), "quadrupole")


def saddle_pair(radius, z_lo, z_hi, opening=2*np.pi/3, arc_segments=6):
    """Two saddles centred on +y and -y (axial legs at 90 deg +/- opening/2 and 270 deg +/- opening/2) with
    axial current density ~ -cos(phi): a positive current gives +y field at the centre."""
    h = opening/2
    paths = []
    for centre in (np.pi/2, 3*np.pi/2):
        lead, trail = centre+h, centre-h           # up the lead leg, across the top, down the trail leg
        if centre > np.pi:
            lead, trail = centre-h, centre+h
        up = np.array([[radius*np.cos(lead), radius*np.sin(lead), z_lo], [radius*np.cos(lead), radius*np.sin(lead), z_hi]])
        top = arc(radius, lead, trail, z_hi, arc_segments)[1:]
        down = np.array([[radius*np.cos(trail), radius*np.sin(trail), z_lo]])
        bottom = arc(radius, trail, lead, z_lo, arc_segments)[1:]
        paths.append(np.vstack((up, top, down, bottom)))
    return paths


def rmf_antenna(frequency, amplitude, parity="even", radius=0.16, half_length=0.6, direction=1, **saddle):
    """Two orthogonal saddle pairs driven I0 cos(wt) (field along x) and I0 sin(wt) (along y): a transverse
    field rotating in +theta for direction=+1. Even parity: one pair per axis over |z| <= half_length.
    Odd parity: each pair split at z = 0 with opposite currents in the two halves, so B_perp reverses sign
    across the midplane. I0 is set so |B_perp| = amplitude at the centre (even) or at z = half_length/2
    (odd). Returns (list of Filaments, I0)."""
    w = 2*np.pi*frequency
    if parity == "even":
        y_pair = saddle_pair(radius, -half_length, half_length, **saddle)
        probe = np.zeros((1, 3))
    elif parity == "odd":
        y_pair = saddle_pair(radius, 0., half_length, **saddle)+[p[::-1] for p in saddle_pair(radius, -half_length, 0., **saddle)]
        probe = np.array([[0., 0., half_length/2]])
    else:
        raise ValueError("parity is 'even' or 'odd'")
    x_pair = rotate_z(y_pair, -np.pi/2)
    b = Filaments(y_pair, None).per_ampere(probe, potential=False)[0][0]
    i0 = amplitude/abs(b[1])*np.sign(b[1])
    pair_x = Filaments(x_pair, lambda t: (i0*np.cos(w*t), -i0*w*np.sin(w*t)), "RMF x")
    pair_y = Filaments(y_pair, lambda t: (direction*i0*np.sin(w*t), direction*i0*w*np.cos(w*t)), "RMF y")
    return [pair_x, pair_y], abs(i0)


# ---------------------------------------------------------------- integrator

@dataclass
class Walls:
    """Absorbing surfaces: r = r_outer, r = r_inner (if > 0) and |z| = z_end."""
    r_outer: float
    r_inner: float = 0.
    z_end: float = np.inf

    def outside(self, X):
        r2 = X[:, 0]**2+X[:, 1]**2
        return (r2 >= self.r_outer**2) | (r2 <= self.r_inner**2) | (np.abs(X[:, 2]) >= self.z_end)

    def crossing(self, X0, X1):
        """(fraction s of the step, surface code) of the first surface crossed on X0 -> X1."""
        d = X1-X0
        n = len(X0)
        best, code = np.full(n, np.inf), np.zeros(n, int)
        a = d[:, 0]**2+d[:, 1]**2
        b = X0[:, 0]*d[:, 0]+X0[:, 1]*d[:, 1]
        for radius, sign, surface in ((self.r_outer, 1, 1), (self.r_inner, -1, 2)):
            if radius <= 0:
                continue
            c = X0[:, 0]**2+X0[:, 1]**2-radius**2
            disc = np.maximum(b*b-a*c, 0.)
            with np.errstate(divide="ignore", invalid="ignore"):
                s = (-b+sign*np.sqrt(disc))/a
            hit = (a > 0) & (s >= 0) & (s <= 1) & (sign*(X1[:, 0]**2+X1[:, 1]**2-radius**2) >= 0)
            take = hit & (s < best)
            best[take], code[take] = s[take], surface
        for z, surface in ((self.z_end, 3), (-self.z_end, 4)):
            if not np.isfinite(z):
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                s = (z-X0[:, 2])/d[:, 2]
            hit = (s >= 0) & (s <= 1) & (np.sign(z)*(X1[:, 2]-z) >= 0)
            take = hit & (s < best)
            best[take], code[take] = s[take], surface
        unresolved = code == 0                       # numerically on a surface: take the step end
        best[unresolved] = 1.
        return best, code


def gamma_of(u):
    return np.sqrt(1+np.einsum("nk,nk->n", u, u)/C0**2)


def boris_kick(u, E, B, qm, dt):
    """Relativistic Boris: half electric kick, magnetic rotation with gamma at the half-kicked momentum,
    half electric kick. u = gamma v."""
    h = 0.5*qm*dt
    um = u+h*E
    t = (h/gamma_of(um))[:, None]*B
    s = 2*t/(1+np.einsum("nk,nk->n", t, t))[:, None]
    up = um+np.cross(um+np.cross(um, t), s)
    return up+h*E


@dataclass
class Population:
    """Particles of one species; positions x [m] and momenta per mass u = gamma v [m/s]."""
    species: str
    x: np.ndarray
    u: np.ndarray
    ids: np.ndarray = None
    lost_time: np.ndarray = None
    lost_position: np.ndarray = None
    lost_surface: np.ndarray = None
    lost_energy_ev: np.ndarray = None
    initial: dict = field(default_factory=dict)

    def __post_init__(self):
        n = len(self.x)
        self.x, self.u = np.array(self.x, float), np.array(self.u, float)
        self.ids = np.arange(n) if self.ids is None else np.asarray(self.ids)
        self.count = n
        self.lost_time = np.full(n, np.nan)
        self.lost_position = np.full((n, 3), np.nan)
        self.lost_surface = np.zeros(n, int)
        self.lost_energy_ev = np.full(n, np.nan)
        self.initial = dict(x=self.x.copy(), u=self.u.copy())

    @property
    def mass(self):
        return SPECIES[self.species][0]

    @property
    def charge(self):
        return SPECIES[self.species][1]

    def energy_ev(self, u=None):
        """Kinetic energy (gamma - 1) m c^2 in eV, as m u^2/(gamma + 1) to keep the digits at low energy."""
        u = self.u if u is None else u
        return self.mass*np.einsum("nk,nk->n", u, u)/(gamma_of(u)+1)/QE


def track(population, fields, t0, dt, steps, walls=None, observe=None):
    """Advance `population` in place by `steps` Boris steps of dt from t0. Lost particles are removed from
    the active arrays (x, u, ids) and recorded by id. observe(t, pop, B_half) is called after every step
    with the active particles and the field B at their half-step positions."""
    p = population
    qm = p.charge/p.mass
    for n in range(steps):
        t = t0+n*dt
        x_half = p.x+0.5*dt*p.u/gamma_of(p.u)[:, None]
        B, E = fields.fields(x_half, t+0.5*dt)
        p.u = boris_kick(p.u, E, B, qm, dt)
        x_new = x_half+0.5*dt*p.u/gamma_of(p.u)[:, None]
        if walls is not None:
            out = walls.outside(x_new)
            if out.any():
                s, code = walls.crossing(p.x[out], x_new[out])
                ids = p.ids[out]
                p.lost_time[ids] = t+s*dt
                p.lost_position[ids] = p.x[out]+s[:, None]*(x_new[out]-p.x[out])
                p.lost_surface[ids] = code
                p.lost_energy_ev[ids] = p.energy_ev(p.u[out])
                keep = ~out
                x_new, B = x_new[keep], B[keep]
                p.u, p.ids = p.u[keep], p.ids[keep]
        p.x = x_new
        if observe is not None:
            observe(t+dt, p, B)
        if not len(p.ids):
            break
    return p


def canonical_angular_momentum(population, fields, t):
    """p_theta per particle: gamma m (x v_y - y v_x) + q psi(x, t) [kg m^2/s]."""
    p = population
    return p.mass*(p.x[:, 0]*p.u[:, 1]-p.x[:, 1]*p.u[:, 0])+p.charge*fields.psi(p.x, t)


def maxwellian(rng, n, temperature_ev, mass):
    """Non-relativistic Maxwellian velocities (n, 3) for temperatures (scalar or per particle) in eV."""
    sigma = np.sqrt(np.asarray(temperature_ev, float)*QE/mass)
    return rng.normal(size=(n, 3))*np.reshape(sigma, (-1, 1) if np.ndim(sigma) else ())


def momentum_from_velocity(v):
    """u = gamma v for velocities v (n, 3) below c."""
    return v/np.sqrt(1-np.einsum("nk,nk->n", v, v)/C0**2)[:, None]
