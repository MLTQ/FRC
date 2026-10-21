"""3D hybrid model in the 12-inch device (M2 of docs/HYBRID_3D.md): gas tube, vacuum gap, conducting coil.

Cells (by centre) are gas (inside the tube, where ions live), vacuum gap (between the tube and the coil, and
beyond the tube ends up to the conducting end planes) or conductor (the coil and end planes, and everything
outside them).
- Faces touching a gas cell are advanced by Faraday's law with the Ohm's-law E on their edges (edges touching
  a gas cell). Edges elsewhere carry no E.
- Faces inside the gap hold the exact vacuum field: J = 0 there, so B = -grad(chi), with chi solving
  div B = 0 in every gap cell (a Neumann problem, data from the gas-side and conductor faces). This is
  re-solved after every field stage. The gap is simply connected (it closes over the tube ends), so no
  harmonic field with circulation is needed. An axial conductor through the device (a rod) would change that.
- Conductor faces keep their normal B (perfect conductor, flux conserver), as in the 2D formation model,
  whose fixed induced-flux boundary sits at the coil.
- Ions that enter a non-gas cell are absorbed (the tube wall and end walls, as a staircase on the cells).
- The Ohm's-law kernels are numba versions of plasma_bench/hybrid3d.py's numpy operators. They are checked
  against them in examples/hybrid3d_verification.py.
"""
from dataclasses import dataclass
import numpy as np
from numba import njit, prange, get_num_threads
from scipy.constants import mu_0 as MU0, e as QE
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu
from .hybrid3d import Box, deposit, boris, HybridParams, divergence_faces

GAS, GAP, CONDUCTOR = 0, 1, 2


@dataclass(frozen=True)
class Device:
    wall_radius: float = .1524        # glass tube inner radius [m] (12" ID)
    tube_half_length: float = 1.1
    coil_radius: float = .165         # single-turn coil (perfect conductor)
    end_gap: float = .02              # vacuum beyond each tube end before the conducting end plane [m]


def device_box(device, dx, dz, periodic_length=None):
    """Cartesian box with at least two conductor layers outside the coil and the end planes. With
    `periodic_length` the tube is infinitely long: the box is that long in z and periodic (no end walls)."""
    half = device.coil_radius+2.5*dx
    nxy = int(np.ceil(2*half/dx))
    if periodic_length is not None:
        nz = int(round(periodic_length/dz))
    else:
        zc = device.tube_half_length+device.end_gap
        nz = int(np.ceil(2*(zc+2.5*dz)/dz))
    hx, hz = nxy*dx/2, nz*dz/2
    return Box((nxy, nxy, nz), (-hx, -hx, -hz), (hx, hx, hz))


def periodic(device):
    """The same cross-section with an infinitely long tube (use with device_box(..., periodic_length=L))."""
    return Device(device.wall_radius, np.inf, device.coil_radius, 0.)


def cell_types(box, device):
    xs = [box.lo[a]+box.d[a]*(np.arange(box.shape[a])+.5) for a in range(3)]
    X, Y, Z = np.meshgrid(*xs, indexing="ij")
    r = np.hypot(X, Y)
    t = np.full(box.shape, GAP, np.int8)
    t[(r < device.wall_radius) & (np.abs(Z) < device.tube_half_length)] = GAS
    # the end planes sit at least one cell beyond the tube ends, so the gap closes over the ends (simply connected)
    t[(r >= device.coil_radius) | (np.abs(Z) > device.tube_half_length+max(device.end_gap, 1.01*box.d[2]))] = CONDUCTOR
    return t


def _shift(a, axis):
    """a at the neighbour on the low side along axis (periodic index; boundaries are conductor)."""
    return np.roll(a, 1, axis=axis)


def masks(types):
    """Face masks (active: touches gas; gap: between two gap cells, set by the projector) and edge masks
    (touch gas, not on a conductor surface, where tangential E = 0). Faces between gap and conductor are frozen
    conductor data."""
    gas = types == GAS
    gap = types == GAP
    face_active, face_gap = [], []
    for a in range(3):
        g = gas | _shift(gas, a)
        face_active.append(g)
        face_gap.append(gap & _shift(gap, a))
    cond = types == CONDUCTOR
    edge = []
    for a in range(3):
        p, q = [b for b in range(3) if b != a]
        g = gas | _shift(gas, p) | _shift(gas, q) | _shift(_shift(gas, p), q)
        c = cond | _shift(cond, p) | _shift(cond, q) | _shift(_shift(cond, p), q)
        edge.append(g & ~c)          # tangential E vanishes on a conductor surface
    return face_active, face_gap, edge


class GapProjector:
    """Exact vacuum field in the gap cells: B = -grad(chi), div B = 0 cell by cell (see module doc)."""
    def __init__(self, types, d):
        self.d = np.asarray(d, float)
        shape = types.shape
        gap = types == GAP
        idx = -np.ones(shape, np.int64)
        idx[gap] = np.arange(int(gap.sum()))
        self.idx, self.n = idx, int(gap.sum())
        if self.n == 0:
            self.known, self.gapfaces, self.pins, self.components = [], [], np.zeros(0, np.int64), 0
            return
        rows, cols, vals = [], [], []
        # known-face contributions to cell rows: rhs[c] -= sign * B_face / d_a
        self.known = []          # (axis, face flat index, gap cell index, sign/d)
        self.gapfaces = []       # (axis, face flat index, low cell index, high cell index)
        flat = np.arange(np.prod(shape)).reshape(shape)
        for a in range(3):
            lo_idx = _shift(idx, a)                # cell on the low side of face [i] along a
            hi_idx = idx                           # cell on the high side
            both = (lo_idx >= 0) & (hi_idx >= 0)
            f = flat[both]
            self.gapfaces.append((a, f, lo_idx[both], hi_idx[both]))
            l, h = lo_idx[both], hi_idx[both]
            w = 1/self.d[a]**2
            rows += [l, h, l, h]; cols += [l, h, h, l]; vals += [np.full(l.size, w)]*2+[np.full(l.size, -w)]*2
            # faces with exactly one gap neighbour carry known B
            only_hi = (hi_idx >= 0) & (lo_idx < 0)     # the face is the gap cell's low face: div gets -B/d
            only_lo = (lo_idx >= 0) & (hi_idx < 0)     # the face is the gap cell's high face: div gets +B/d
            self.known.append((a, flat[only_hi], hi_idx[only_hi], np.full(only_hi.sum(), -1/self.d[a])))
            self.known.append((a, flat[only_lo], lo_idx[only_lo], np.full(only_lo.sum(), 1/self.d[a])))
        A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(self.n, self.n)).tocsr()
        ncomp, labels = connected_components(A, directed=False)
        pins = [int(np.nonzero(labels == c)[0][0]) for c in range(ncomp)]
        A = A.tolil()
        for p in pins:
            A[p, :] = 0
            A[p, p] = 1.
        self.pins = np.array(pins, np.int64)
        self.components = ncomp
        self.lu = splu(A.tocsc())

    def __call__(self, b):
        """Overwrite the gap faces of b (list of three face arrays) with the vacuum field."""
        if self.n == 0:
            return b
        rhs = np.zeros(self.n)
        for a, f, c, s in self.known:
            np.add.at(rhs, c, -s*b[a].ravel()[f])
        rhs[self.pins] = 0.
        chi = self.lu.solve(rhs)
        for a, f, lo, hi in self.gapfaces:
            b[a].ravel()[f] = -(chi[hi]-chi[lo])/self.d[a]
        return b


@njit(parallel=True, cache=True)
def curl_faces_nb(bx, by, bz, dx, dy, dz, jx, jy, jz):
    nx, ny, nz = bx.shape
    for i in prange(nx):
        im = (i-1) % nx
        for j in range(ny):
            jm = (j-1) % ny
            for k in range(nz):
                km = (k-1) % nz
                jx[i, j, k] = ((bz[i, j, k]-bz[i, jm, k])/dy-(by[i, j, k]-by[i, j, km])/dz)/MU0
                jy[i, j, k] = ((bx[i, j, k]-bx[i, j, km])/dz-(bz[i, j, k]-bz[im, j, k])/dx)/MU0
                jz[i, j, k] = ((by[i, j, k]-by[im, j, k])/dx-(bx[i, j, k]-bx[i, jm, k])/dy)/MU0


@njit(parallel=True, cache=True)
def ohm_nodes_nb(bx, by, bz, jx, jy, jz, n, fx, fy, fz, pe, floor, q, dx, dy, dz, en, bn, jn):
    """Node E (ideal + Hall + pressure), node B and node J."""
    nx, ny, nz = bx.shape
    for i in prange(nx):
        im = (i-1) % nx; ip = (i+1) % nx
        for j in range(ny):
            jm = (j-1) % ny; jp = (j+1) % ny
            for k in range(nz):
                km = (k-1) % nz; kp = (k+1) % nz
                bnx = .25*(bx[i, j, k]+bx[i, jm, k]+bx[i, j, km]+bx[i, jm, km])
                bny = .25*(by[i, j, k]+by[im, j, k]+by[i, j, km]+by[im, j, km])
                bnz = .25*(bz[i, j, k]+bz[im, j, k]+bz[i, jm, k]+bz[im, jm, k])
                jnx = .5*(jx[i, j, k]+jx[im, j, k]); jny = .5*(jy[i, j, k]+jy[i, jm, k]); jnz = .5*(jz[i, j, k]+jz[i, j, km])
                ne = max(n[i, j, k], floor)
                inv = 1./(ne*q)
                ux = fx[i, j, k]/ne; uy = fy[i, j, k]/ne; uz = fz[i, j, k]/ne
                ex = -(uy*bnz-uz*bny)+(jny*bnz-jnz*bny)*inv
                ey = -(uz*bnx-ux*bnz)+(jnz*bnx-jnx*bnz)*inv
                ez = -(ux*bny-uy*bnx)+(jnx*bny-jny*bnx)*inv
                ex -= (pe[ip, j, k]-pe[im, j, k])/(2*dx)*inv
                ey -= (pe[i, jp, k]-pe[i, jm, k])/(2*dy)*inv
                ez -= (pe[i, j, kp]-pe[i, j, km])/(2*dz)*inv
                en[0, i, j, k] = ex; en[1, i, j, k] = ey; en[2, i, j, k] = ez
                bn[0, i, j, k] = bnx; bn[1, i, j, k] = bny; bn[2, i, j, k] = bnz
                jn[0, i, j, k] = jnx; jn[1, i, j, k] = jny; jn[2, i, j, k] = jnz


@njit(parallel=True, cache=True)
def ohm_edges_nb(en, jx, jy, jz, eta, eta_h, dx, dy, dz, mx, my, mz, ex, ey, ez):
    """Edge E = average of node E along the edge + eta J - eta_H lap J, masked to edges touching gas."""
    nx, ny, nz = jx.shape
    for i in prange(nx):
        im = (i-1) % nx; ip = (i+1) % nx
        for j in range(ny):
            jm = (j-1) % ny; jp = (j+1) % ny
            for k in range(nz):
                km = (k-1) % nz; kp = (k+1) % nz
                vx = vy = vz = 0.
                if mx[i, j, k]:
                    vx = .5*(en[0, i, j, k]+en[0, ip, j, k])+eta*jx[i, j, k]
                    if eta_h > 0:
                        vx -= eta_h*((jx[ip, j, k]-2*jx[i, j, k]+jx[im, j, k])/dx**2+(jx[i, jp, k]-2*jx[i, j, k]+jx[i, jm, k])/dy**2
                                     + (jx[i, j, kp]-2*jx[i, j, k]+jx[i, j, km])/dz**2)
                if my[i, j, k]:
                    vy = .5*(en[1, i, j, k]+en[1, i, jp, k])+eta*jy[i, j, k]
                    if eta_h > 0:
                        vy -= eta_h*((jy[ip, j, k]-2*jy[i, j, k]+jy[im, j, k])/dx**2+(jy[i, jp, k]-2*jy[i, j, k]+jy[i, jm, k])/dy**2
                                     + (jy[i, j, kp]-2*jy[i, j, k]+jy[i, j, km])/dz**2)
                if mz[i, j, k]:
                    vz = .5*(en[2, i, j, k]+en[2, i, j, kp])+eta*jz[i, j, k]
                    if eta_h > 0:
                        vz -= eta_h*((jz[ip, j, k]-2*jz[i, j, k]+jz[im, j, k])/dx**2+(jz[i, jp, k]-2*jz[i, j, k]+jz[i, jm, k])/dy**2
                                     + (jz[i, j, kp]-2*jz[i, j, k]+jz[i, j, km])/dz**2)
                ex[i, j, k] = vx; ey[i, j, k] = vy; ez[i, j, k] = vz


@njit(parallel=True, cache=True)
def faraday_nb(ex, ey, ez, dx, dy, dz, ax, ay, az, rx, ry, rz):
    """dB/dt = -curl E on active faces (zero elsewhere)."""
    nx, ny, nz = ex.shape
    for i in prange(nx):
        ip = (i+1) % nx
        for j in range(ny):
            jp = (j+1) % ny
            for k in range(nz):
                kp = (k+1) % nz
                rx[i, j, k] = -((ez[i, jp, k]-ez[i, j, k])/dy-(ey[i, j, kp]-ey[i, j, k])/dz) if ax[i, j, k] else 0.
                ry[i, j, k] = -((ex[i, j, kp]-ex[i, j, k])/dz-(ez[ip, j, k]-ez[i, j, k])/dx) if ay[i, j, k] else 0.
                rz[i, j, k] = -((ey[ip, j, k]-ey[i, j, k])/dx-(ex[i, jp, k]-ex[i, j, k])/dy) if az[i, j, k] else 0.


@njit(parallel=True, cache=True)
def inside_gas(x, lo, d, shape, types, alive):
    """alive[p] = particle p sits in a gas cell (staircase tube wall and ends absorb the rest)."""
    for p in prange(x.shape[0]):
        i = int(np.floor((x[p, 0]-lo[0])/d[0])); j = int(np.floor((x[p, 1]-lo[1])/d[1])); k = int(np.floor((x[p, 2]-lo[2])/d[2]))
        ok = 0 <= i < shape[0] and 0 <= j < shape[1] and 0 <= k < shape[2]
        alive[p] = ok and types[i, j, k] == 0


class DeviceHybrid:
    """Kinetic ions + electron fluid in the device. Particles: x at t, v at t - dt/2, weights w."""
    def __init__(self, box, device, params, b_faces, x, v, w, dt, n_ref=None):
        self.box, self.device, self.p, self.dt = box, device, params, float(dt)
        self.d = box.d
        self.lo = np.asarray(box.lo, float); self.len = np.asarray(box.hi, float)-self.lo
        self.shape = np.asarray(box.shape, np.int64)
        self.types = cell_types(box, device)
        self.face_active, self.face_gap, self.edge_mask = masks(self.types)
        self.u8 = [m.astype(np.uint8) for m in self.edge_mask]
        # In a z-periodic tube the gap is a ring around the plasma: besides gradients it carries one harmonic
        # field, a uniform Bz (axial flux through the ring). Its gap z-faces are advanced by Faraday (E = 0 on
        # gap-interior edges, Ohm's law on the plasma surface), so the ring flux changes by the line integral of E
        # around the plasma, and the projection keeps their mean. (A ring circulation B_theta would need a second
        # harmonic; no axial current threads the plasma in these runs.)
        self.periodic_z = not np.isfinite(device.tube_half_length)
        faraday_faces = [m.copy() for m in self.face_active]
        if self.periodic_z:
            faraday_faces[2] |= self.face_gap[2]
        self.f8 = [m.astype(np.uint8) for m in faraday_faces]
        self._gap_project = GapProjector(self.types, self.d)
        self.b = self.project([np.array(c, float) for c in b_faces])
        self.x, self.v = np.ascontiguousarray(x, float), np.ascontiguousarray(v, float)
        self.w = np.broadcast_to(np.asarray(w, float), (self.x.shape[0],)).copy()
        self.qm = params.ion_charge/params.ion_mass
        self.nthreads = get_num_threads()
        self._node_inv = self._node_fraction()
        self._absorb()
        mom = self._moments(self.x)
        self.n = mom[0]
        self.flux_half = mom[1:]
        self.flux_prev = mom[1:].copy()
        self.n_ref = float(self.n[self.types == GAS].mean()) if n_ref is None else n_ref
        self.time, self.steps, self.substeps = 0., 0, 0
        self.lost = dict(count=0., energy_j=0., ends=0., wall=0.)
        shape = tuple(box.shape)
        self._j = [np.zeros(shape) for _ in range(3)]
        self._e = [np.zeros(shape) for _ in range(3)]
        self._en = np.zeros((3,)+shape); self._bn = np.zeros((3,)+shape); self._jn = np.zeros((3,)+shape)
        self._rate = [np.zeros(shape) for _ in range(3)]

    def project(self, b):
        if self.periodic_z and self.face_gap[2].any():
            ring = float(b[2][self.face_gap[2]].mean())
            b = self._gap_project(b)
            b[2][self.face_gap[2]] += ring-float(b[2][self.face_gap[2]].mean())
            return b
        return self._gap_project(b)

    def _moments(self, x, v=None):
        """Node density and flux, divided by the fraction of each node's 8 surrounding cells that are gas, so
        a uniform plasma has uniform node density up to the staircase wall (as the 2D model's wall nodes carry
        half dual areas). Nodes with no gas cell get zero."""
        mom = deposit(np.ascontiguousarray(x), self.v if v is None else v, self.w, self.lo, self.d, self.shape, self.nthreads)
        return mom*self._node_inv

    def _node_fraction(self):
        gas = (self.types == GAS).astype(float)
        f = gas.copy()
        for a in range(3):
            f = .5*(f+np.roll(f, 1, axis=a))      # node (i,j,k) averages cells (i-1..i, j-1..j, k-1..k)
        inv = np.zeros_like(f)
        inv[f > 0] = 1/f[f > 0]
        self.node_gas_fraction = f
        return inv

    def _absorb(self):
        alive = np.zeros(self.x.shape[0], np.bool_)
        inside_gas(self.x, self.lo, self.d, self.shape, self.types, alive)
        if not alive.all():
            dead = ~alive
            self.lost["count"] += float(self.w[dead].sum())
            self.lost["energy_j"] += float(.5*self.p.ion_mass*np.sum(self.w[dead]*np.sum(self.v[dead]**2, axis=1)))
            ends = np.abs(self.x[dead, 2]) >= self.device.tube_half_length
            self.lost["ends"] += float(self.w[dead][ends].sum()); self.lost["wall"] += float(self.w[dead][~ends].sum())
            self.x, self.v, self.w = self.x[alive], self.v[alive], self.w[alive]

    def _pressure(self, n):
        te = self.p.te_ev*QE
        if te == 0:
            return np.zeros_like(n)
        if self.p.gamma_e is None:
            return n*te
        return n*te*(np.maximum(n, 1e-300)/self.n_ref)**(self.p.gamma_e-1)

    def electric(self, b, n, flux):
        """Fill edge E (masked) and node E/B/J."""
        curl_faces_nb(b[0], b[1], b[2], *self.d, *self._j)
        ohm_nodes_nb(b[0], b[1], b[2], *self._j, n, flux[0], flux[1], flux[2], self._pressure(n), self.p.n_floor,
                     self.p.ion_charge, *self.d, self._en, self._bn, self._jn)
        ohm_edges_nb(self._en, *self._j, self.p.eta, self.p.eta_hyper, *self.d, *self.u8, *self._e)
        return self._e

    def rate(self, b, n, flux):
        e = self.electric(b, n, flux)
        faraday_nb(e[0], e[1], e[2], *self.d, *self.f8, *self._rate)
        return [r.copy() for r in self._rate]

    def field_substeps(self):
        live = self.node_gas_fraction > 0
        b = max(float(np.sqrt(sum(np.max(np.abs(c[self.face_active[a]]))**2 for a, c in enumerate(self.b)))), 1e-300)
        nmin = max(float(self.n[live].min()) if live.any() else 0., self.p.n_floor, 1e-300)
        kmax = np.pi*np.sqrt(np.sum(1/self.d**2))
        va = b/np.sqrt(MU0*nmin*self.p.ion_mass)
        omega = kmax*va+kmax**2*b/(MU0*nmin*self.p.ion_charge)+(self.p.eta/MU0)*kmax**2+(self.p.eta_hyper/MU0)*kmax**4
        if self.p.te_ev > 0:
            omega += kmax*np.sqrt(self.p.te_ev*QE*(self.p.gamma_e or 1.)/self.p.ion_mass)
        return max(1, int(np.ceil(self.dt*omega/(self.p.field_cfl*np.sqrt(3)))))

    def step(self):
        dt = self.dt
        flux_now = 1.5*self.flux_half-.5*self.flux_prev
        self.electric(self.b, self.n, flux_now)
        e_nodes = self._en+self.p.eta*self._jn if self.p.eta > 0 else self._en
        x_old = self.x.copy()
        boris(self.x, self.v, np.ascontiguousarray(e_nodes), np.ascontiguousarray(self._bn), self.qm, dt,
              self.lo, self.d, self.shape, self.len)
        mid = np.ascontiguousarray(x_old+.5*dt*self.v)
        self._absorb_with(mid)
        n_new = self._moments(self.x)[0]
        flux_half = self._moments(self._mid)[1:]
        self.flux_prev, self.flux_half = self.flux_half, flux_half
        m = self.field_substeps()
        h = dt/m
        n_old = self.n
        for s in range(m):
            f0 = s*h/dt
            b = self.b
            k1 = self.rate(b, n_old+(n_new-n_old)*f0, flux_half)
            b1 = self.project([b[a]+h*k1[a] for a in range(3)])
            k2 = self.rate(b1, n_old+(n_new-n_old)*(f0+h/dt), flux_half)
            b2 = self.project([.75*b[a]+.25*(b1[a]+h*k2[a]) for a in range(3)])
            k3 = self.rate(b2, n_old+(n_new-n_old)*(f0+.5*h/dt), flux_half)
            self.b = self.project([b[a]/3+2/3*(b2[a]+h*k3[a]) for a in range(3)])
        self.n = n_new
        self.time += dt; self.steps += 1; self.substeps += m
        return m

    def _absorb_with(self, mid):
        """Absorb particles that left the gas; keep the matching mid positions for the flux deposit."""
        alive = np.zeros(self.x.shape[0], np.bool_)
        inside_gas(self.x, self.lo, self.d, self.shape, self.types, alive)
        if not alive.all():
            dead = ~alive
            self.lost["count"] += float(self.w[dead].sum())
            self.lost["energy_j"] += float(.5*self.p.ion_mass*np.sum(self.w[dead]*np.sum(self.v[dead]**2, axis=1)))
            ends = np.abs(self.x[dead, 2]) >= self.device.tube_half_length
            self.lost["ends"] += float(self.w[dead][ends].sum()); self.lost["wall"] += float(self.w[dead][~ends].sum())
            self.x, self.v, self.w, mid = self.x[alive], self.v[alive], self.w[alive], mid[alive]
        self._mid = np.ascontiguousarray(mid)

    # diagnostics
    def magnetic_energy(self):
        return float(sum(np.sum(c**2) for c in self.b)/(2*MU0)*self.box.cell_volume)

    def kinetic_energy(self):
        return float(.5*self.p.ion_mass*np.sum(self.w*np.sum(self.v**2, axis=1)))

    def particles(self):
        return float(self.w.sum())

    def max_divergence(self, region=None):
        div = divergence_faces(*self.b, self.d)
        sel = self.types != CONDUCTOR if region is None else region
        scale = max(float(np.sqrt(sum(np.max(np.abs(c))**2 for c in self.b))), 1e-300)/self.d.min()
        return float(np.abs(div[sel]).max()/scale)
