"""3D hybrid plasma model: kinetic full-orbit ions, massless quasi-neutral electron fluid, Darwin fields.

Yee grid on a periodic Cartesian box (M1 of docs/HYBRID_3D.md; device geometry follows in M2):
- B on faces: Bx[i,j,k] at (i, j+1/2, k+1/2), By at (i+1/2, j, k+1/2), Bz at (i+1/2, j+1/2, k) [cell units];
- E and J on edges: Ex[i,j,k] at (i+1/2, j, k), Ey at (i, j+1/2, k), Ez at (i, j, k+1/2);
- nodes (i, j, k) carry the ion moments and the Ohm's-law field used by the particles.

Electron Ohm's law (massless, n_e = n_i, Darwin J = curl B / mu0):
    E = -u_i x B + J x B/(n e) - grad p_e/(n e) + eta J - eta_H lap J
The ideal, Hall and pressure terms are evaluated at nodes and averaged to edges; eta J and the declared
hyper-resistivity eta_H act on the edge current directly. Faraday dB/dt = -curl E keeps div B at round-off.
"""
from dataclasses import dataclass
import numpy as np
from numba import njit, prange, get_num_threads
from scipy.constants import mu_0 as MU0, e as QE

AXES = (0, 1, 2)


@dataclass(frozen=True)
class Box:
    """Periodic box: `shape` cells between `lo` and `hi` [m]."""
    shape: tuple
    lo: tuple
    hi: tuple

    @property
    def d(self):
        return (np.asarray(self.hi, float)-np.asarray(self.lo, float))/np.asarray(self.shape)

    @property
    def cell_volume(self):
        return float(np.prod(self.d))

    def nodes(self, axis):
        return self.lo[axis]+self.d[axis]*np.arange(self.shape[axis])


def _r(a, shift, axis):
    return np.roll(a, shift, axis=axis)


def curl_edges(ex, ey, ez, d):
    """curl of an edge field, on faces."""
    dx, dy, dz = d
    return ((_r(ez, -1, 1)-ez)/dy-(_r(ey, -1, 2)-ey)/dz,
            (_r(ex, -1, 2)-ex)/dz-(_r(ez, -1, 0)-ez)/dx,
            (_r(ey, -1, 0)-ey)/dx-(_r(ex, -1, 1)-ex)/dy)


def curl_faces(bx, by, bz, d):
    """curl of a face field, on edges."""
    dx, dy, dz = d
    return ((bz-_r(bz, 1, 1))/dy-(by-_r(by, 1, 2))/dz,
            (bx-_r(bx, 1, 2))/dz-(bz-_r(bz, 1, 0))/dx,
            (by-_r(by, 1, 0))/dx-(bx-_r(bx, 1, 1))/dy)


def divergence_faces(bx, by, bz, d):
    """div of a face field, at cell centres."""
    return (_r(bx, -1, 0)-bx)/d[0]+(_r(by, -1, 1)-by)/d[1]+(_r(bz, -1, 2)-bz)/d[2]


def faces_to_nodes(bx, by, bz):
    def avg(a, p, q):
        return .25*(a+_r(a, 1, p)+_r(a, 1, q)+_r(_r(a, 1, p), 1, q))
    return np.stack((avg(bx, 1, 2), avg(by, 0, 2), avg(bz, 0, 1)))


def edges_to_nodes(ex, ey, ez):
    return np.stack((.5*(ex+_r(ex, 1, 0)), .5*(ey+_r(ey, 1, 1)), .5*(ez+_r(ez, 1, 2))))


def nodes_to_edges(v):
    return tuple(.5*(v[a]+_r(v[a], -1, a)) for a in AXES)


def laplacian(f, d):
    return sum((_r(f, -1, a)-2*f+_r(f, 1, a))/d[a]**2 for a in AXES)


def gradient_nodes(p, d):
    return np.stack([(_r(p, -1, a)-_r(p, 1, a))/(2*d[a]) for a in AXES])


@njit(parallel=True, cache=True)
def deposit(x, v, w, lo, d, shape, nthreads):
    """Trilinear (CIC) node deposit of number density and particle flux n u [1/m^3, 1/(m^2 s)]."""
    nx, ny, nz = shape[0], shape[1], shape[2]
    local = np.zeros((nthreads, 4, nx, ny, nz))
    n = x.shape[0]
    chunk = (n+nthreads-1)//nthreads
    for t in prange(nthreads):
        for p in range(t*chunk, min(n, (t+1)*chunk)):
            sx = (x[p, 0]-lo[0])/d[0]; sy = (x[p, 1]-lo[1])/d[1]; sz = (x[p, 2]-lo[2])/d[2]
            fx = np.floor(sx); fy = np.floor(sy); fz = np.floor(sz)
            i0 = int(fx) % nx; j0 = int(fy) % ny; k0 = int(fz) % nz
            ax = sx-fx; ay = sy-fy; az = sz-fz
            i1 = (i0+1) % nx; j1 = (j0+1) % ny; k1 = (k0+1) % nz
            vx = v[p, 0]; vy = v[p, 1]; vz = v[p, 2]; wp = w[p]
            for c in range(8):
                ii = i1 if c & 1 else i0
                jj = j1 if (c >> 1) & 1 else j0
                kk = k1 if (c >> 2) & 1 else k0
                wt = wp*(ax if c & 1 else 1-ax)*(ay if (c >> 1) & 1 else 1-ay)*(az if (c >> 2) & 1 else 1-az)
                local[t, 0, ii, jj, kk] += wt
                local[t, 1, ii, jj, kk] += wt*vx
                local[t, 2, ii, jj, kk] += wt*vy
                local[t, 3, ii, jj, kk] += wt*vz
    out = np.zeros((4, nx, ny, nz))
    for t in range(nthreads):
        out += local[t]
    return out/(d[0]*d[1]*d[2])


@njit(parallel=True, cache=True)
def boris(x, v, e_nodes, b_nodes, qm, dt, lo, d, shape, length):
    """Boris push with node fields gathered trilinearly at x; positions wrap periodically."""
    nx, ny, nz = shape[0], shape[1], shape[2]
    h = .5*qm*dt
    for p in prange(x.shape[0]):
        sx = (x[p, 0]-lo[0])/d[0]; sy = (x[p, 1]-lo[1])/d[1]; sz = (x[p, 2]-lo[2])/d[2]
        fx = np.floor(sx); fy = np.floor(sy); fz = np.floor(sz)
        i0 = int(fx) % nx; j0 = int(fy) % ny; k0 = int(fz) % nz
        ax = sx-fx; ay = sy-fy; az = sz-fz
        i1 = (i0+1) % nx; j1 = (j0+1) % ny; k1 = (k0+1) % nz
        ex = ey = ez = bx = by = bz = 0.
        for c in range(8):
            ii = i1 if c & 1 else i0
            jj = j1 if (c >> 1) & 1 else j0
            kk = k1 if (c >> 2) & 1 else k0
            wt = (ax if c & 1 else 1-ax)*(ay if (c >> 1) & 1 else 1-ay)*(az if (c >> 2) & 1 else 1-az)
            ex += wt*e_nodes[0, ii, jj, kk]; ey += wt*e_nodes[1, ii, jj, kk]; ez += wt*e_nodes[2, ii, jj, kk]
            bx += wt*b_nodes[0, ii, jj, kk]; by += wt*b_nodes[1, ii, jj, kk]; bz += wt*b_nodes[2, ii, jj, kk]
        ux = v[p, 0]+h*ex; uy = v[p, 1]+h*ey; uz = v[p, 2]+h*ez
        tx = h*bx; ty = h*by; tz = h*bz
        f = 2./(1.+tx*tx+ty*ty+tz*tz)
        sx2 = f*tx; sy2 = f*ty; sz2 = f*tz
        px = ux+(uy*tz-uz*ty); py = uy+(uz*tx-ux*tz); pz = uz+(ux*ty-uy*tx)
        ux += py*sz2-pz*sy2; uy += pz*sx2-px*sz2; uz += px*sy2-py*sx2
        v[p, 0] = ux+h*ex; v[p, 1] = uy+h*ey; v[p, 2] = uz+h*ez
        x[p, 0] = lo[0]+(x[p, 0]+dt*v[p, 0]-lo[0]) % length[0]
        x[p, 1] = lo[1]+(x[p, 1]+dt*v[p, 1]-lo[1]) % length[1]
        x[p, 2] = lo[2]+(x[p, 2]+dt*v[p, 2]-lo[2]) % length[2]


@dataclass(frozen=True)
class HybridParams:
    ion_mass: float
    ion_charge: float = QE
    te_ev: float = 0.             # electron temperature [eV]; 0 = cold electrons
    gamma_e: float | None = None  # None: isothermal electrons; else adiabatic p_e = n T_e (n/n0)^(gamma-1)
    eta: float = 0.               # physical resistivity [ohm m]
    eta_hyper: float = 0.         # DECLARED hyper-resistivity [ohm m^3]: damps grid-scale modes; its heat is reported
    n_floor: float = 0.           # DECLARED density floor in Ohm's law [1/m^3]
    field_cfl: float = .5         # fraction of the SSP-RK3 imaginary-axis limit sqrt(3) for the field subcycles


class Hybrid:
    """Kinetic ions + electron fluid on a periodic Yee box. Particles: x at t, v at t - dt/2, weight w
    (physical ions per macro-particle)."""
    def __init__(self, box, params, b_faces, x, v, w, dt, n_ref=None):
        self.box, self.p, self.dt = box, params, float(dt)
        self.d = box.d
        self.lo = np.asarray(box.lo, float); self.len = np.asarray(box.hi, float)-self.lo
        self.shape = np.asarray(box.shape, np.int64)
        self.b = [np.array(c, float) for c in b_faces]
        self.x, self.v = np.ascontiguousarray(x, float), np.ascontiguousarray(v, float)
        self.w = np.broadcast_to(np.asarray(w, float), (self.x.shape[0],)).copy()
        self.qm = params.ion_charge/params.ion_mass
        self.nthreads = get_num_threads()
        mom = deposit(self.x, self.v, self.w, self.lo, self.d, self.shape, self.nthreads)
        self.n = mom[0]
        self.flux_half = mom[1:]           # n u at t - dt/2 (approximated at x(t) for the first step)
        self.flux_prev = mom[1:].copy()
        self.n_ref = float(np.mean(self.n)) if n_ref is None else n_ref
        self.time, self.steps, self.substeps = 0., 0, 0
        self.hyper_heat = self.ohmic_heat = 0.

    # --- Ohm's law -------------------------------------------------------------------------
    def _pressure(self, n):
        te = self.p.te_ev*QE
        if te == 0:
            return None
        if self.p.gamma_e is None:
            return n*te
        return n*te*(np.maximum(n, 1e-300)/self.n_ref)**(self.p.gamma_e-1)

    def electric(self, b, n, flux):
        """Edge E and node E (for the particles) from Ohm's law; returns (edges, nodes, edge J)."""
        neff = np.maximum(n, self.p.n_floor) if self.p.n_floor > 0 else n
        j_edges = [c/MU0 for c in curl_faces(*b, self.d)]
        bn = faces_to_nodes(*b)
        jn = edges_to_nodes(*j_edges)
        u = flux/neff
        ne = neff*self.p.ion_charge
        e = -np.cross(u, bn, axis=0)+np.cross(jn, bn, axis=0)/ne
        pe = self._pressure(n)
        if pe is not None:
            e -= gradient_nodes(pe, self.d)/ne
        edges = list(nodes_to_edges(e))
        if self.p.eta > 0:
            for a in AXES:
                edges[a] = edges[a]+self.p.eta*j_edges[a]
            e = e+self.p.eta*jn
        if self.p.eta_hyper > 0:
            for a in AXES:
                edges[a] = edges[a]-self.p.eta_hyper*laplacian(j_edges[a], self.d)
        return edges, e, j_edges, bn

    def field_substeps(self):
        b = np.sqrt(sum(np.max(np.abs(c))**2 for c in self.b))
        nmin = max(float(self.n.min()), self.p.n_floor, 1e-300)
        kmax = np.pi*np.sqrt(np.sum(1/self.d**2))
        va = b/np.sqrt(MU0*nmin*self.p.ion_mass)
        omega = kmax*va+kmax**2*b/(MU0*nmin*self.p.ion_charge)+(self.p.eta/MU0)*kmax**2+(self.p.eta_hyper/MU0)*kmax**4
        if self.p.te_ev > 0:
            omega += kmax*np.sqrt(self.p.te_ev*QE*(self.p.gamma_e or 1.)/self.p.ion_mass)
        return max(1, int(np.ceil(self.dt*omega/(self.p.field_cfl*np.sqrt(3)))))

    # --- one ion step -------------------------------------------------------------------------
    def step(self):
        dt = self.dt
        flux_now = 1.5*self.flux_half-.5*self.flux_prev
        _, e_nodes, _, b_nodes = self.electric(self.b, self.n, flux_now)
        x_old = self.x.copy()
        boris(self.x, self.v, np.ascontiguousarray(e_nodes), np.ascontiguousarray(b_nodes), self.qm, dt,
              self.lo, self.d, self.shape, self.len)
        n_new = deposit(self.x, self.v, self.w, self.lo, self.d, self.shape, self.nthreads)[0]
        # flux at the half step: positions midway (unwrapped displacement dt v)
        mid = self.lo+(x_old+.5*dt*self.v-self.lo) % self.len
        flux_half = deposit(np.ascontiguousarray(mid), self.v, self.w, self.lo, self.d, self.shape, self.nthreads)[1:]
        self.flux_prev, self.flux_half = self.flux_half, flux_half
        m = self.field_substeps()
        h = dt/m
        n_old = self.n
        for s in range(m):
            t0 = s*h/dt
            self.b = self._rk3(self.b, h, lambda frac: n_old+(n_new-n_old)*frac, t0, h/dt, flux_half)
        self.n = n_new
        self.time += dt; self.steps += 1; self.substeps += m
        return m

    def _rk3(self, b, h, density, t0, fh, flux):
        def rate(bb, frac):
            edges, _, j, _ = self.electric(bb, density(frac), flux)
            return [-c for c in curl_edges(*edges, self.d)]
        k1 = rate(b, t0)
        b1 = [b[a]+h*k1[a] for a in AXES]
        k2 = rate(b1, t0+fh)
        b2 = [.75*b[a]+.25*(b1[a]+h*k2[a]) for a in AXES]
        k3 = rate(b2, t0+.5*fh)
        return [b[a]/3+2/3*(b2[a]+h*k3[a]) for a in AXES]

    # --- diagnostics ----------------------------------------------------------------------------
    def magnetic_energy(self):
        return float(sum(np.sum(c**2) for c in self.b)/(2*MU0)*self.box.cell_volume)

    def kinetic_energy(self, v=None):
        v = self.v if v is None else v
        return float(.5*self.p.ion_mass*np.sum(self.w*np.sum(v**2, axis=1)))

    def max_divergence(self):
        div = divergence_faces(*self.b, self.d)
        scale = max(np.sqrt(sum(np.max(np.abs(c))**2 for c in self.b)), 1e-300)/self.d.min()
        return float(np.abs(div).max()/scale)
