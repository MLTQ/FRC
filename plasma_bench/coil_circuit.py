"""Lumped pulsed-power circuit for single-turn theta-pinch coils: strap inductances, capacitor modules,
switch timing and crowbars. Linear and explicit about its idealizations; see coil_circuit.md."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import ellipk, ellipe
from .electromagnetics import loop_field
from .coils import Coil
from .waveforms import Callback

STRIP_GMD = np.exp(-1.5)  # geometric mean distance of a thin strip of width d from itself, in units of d


@dataclass(frozen=True)
class Strap:
    """Thin conducting band: inner radius, axial centre, axial width [m].

    `section` is the circuit branch (one pair of terminals). Within a section, straps with the same `turn`
    are in parallel between the same collector plates, so each of their filaments links equal flux; the
    turns of a section are in series (same current, flux linkages add). A single-turn sheet coil is one
    section with every strap on turn 0; a wound mirror coil is one section with one strap per turn."""
    radius: float
    z: float
    width: float
    section: int = 0
    turn: int = 0

    def __post_init__(self):
        if not (np.isfinite(self.radius) and self.radius > 0 and np.isfinite(self.z) and self.width > 0):
            raise ValueError("Strap radius and width must be positive, centre finite")
        if not all(isinstance(v, int) and v >= 0 for v in (self.section, self.turn)):
            raise ValueError("Section and turn indices must be nonnegative integers")


def _log_gmd(lo1, hi1, lo2, hi2):
    """ln of the geometric mean distance between collinear segments [lo1,hi1] and [lo2,hi2]:
    (1/(w1 w2)) double integral of ln|x-y|, with G(u) = u^2 ln|u|/2 - 3u^2/4 (G(0) = 0)."""
    def g(u):
        u = np.asarray(u, dtype=float)
        return np.where(u == 0, 0., .5*u*u*np.log(np.abs(np.where(u == 0, 1., u)))-.75*u*u)
    total = g(hi1-lo2)-g(hi1-hi2)-g(lo1-lo2)+g(lo1-hi2)
    return total/((hi1-lo1)*(hi2-lo2))


def filament_inductance(radii, z, widths):
    """Partial inductance matrix [H] of coaxial thin-strip filaments (uniform current within each).

    Off-diagonal: exact mutual inductance of two circular filaments, 2 pi psi. For strips on the same
    radius the logarithmic near field is replaced by its strip average, M += mu0 a ln(d/g) with g the
    geometric mean distance of the two strips. Diagonal: L = mu0 a (ln(8a/g) - 2), g = d exp(-3/2).
    Both corrections are exact to O((width/a)^2); filament-pair fields elsewhere are exact."""
    radii, z, widths = (np.asarray(v, dtype=float) for v in (radii, z, widths))
    n = radii.size
    m = np.empty((n, n))
    for j in range(n):
        others = np.arange(n) != j
        m[others, j] = 2*np.pi*loop_field(radii[others], z[others], radii[j], z[j])[0]
        same = others & (np.abs(radii-radii[j]) <= 1e-12*radii[j])
        if same.any():
            d = np.abs(z[same]-z[j])
            log_g = _log_gmd(z[same]-widths[same]/2, z[same]+widths[same]/2, z[j]-widths[j]/2, z[j]+widths[j]/2)
            m[same, j] += MU0*radii[j]*(np.log(d)-log_g)
    m[np.diag_indices(n)] = MU0*radii*(np.log(8*radii/(widths*STRIP_GMD))-2)
    return .5*(m+m.T)


def coax_inductance(length, inner, outer):
    """Inductance of a coaxial line [H]: L = mu0 l/(2 pi) ln(b/a) for a conductor of radius a returning through a
    coaxial cage of radius b (thin-walled; field mu0 I/(2 pi r) between them, none outside the cage)."""
    if not (length > 0 and 0 < inner < outer):
        raise ValueError("Coaxial line needs length > 0 and 0 < inner < outer radius")
    return MU0*length/(2*np.pi)*np.log(outer/inner)


def nagaoka_inductance(radius, length, turns=1):
    """Exact inductance of a thin cylindrical current sheet with UNIFORM current (Nagaoka coefficient):
    L = mu0 pi a^2 N^2/l * (4/(3 pi k'))[(k'^2/k^2)(K-E) + E - k], k = 2a/sqrt(4a^2+l^2). Independent check."""
    k = 2*radius/np.hypot(2*radius, length)
    kp = length/np.hypot(2*radius, length)
    m = k*k
    big_k, big_e = ellipk(m), ellipe(m)
    coefficient = 4/(3*np.pi*kp)*((kp*kp/m)*(big_k-big_e)+big_e-k)
    return MU0*np.pi*radius**2*turns**2/length*coefficient


@dataclass
class SectionCoil:
    """Coil reduced to sections. Each filament of a section links the section's flux (perfect conductor on
    the circuit time scale). inductance: section matrix L_ss' [H]; distribution: filament current per unit
    section current (n_filaments x n_sections; neighbouring sections redistribute current, totals fixed)."""
    straps: tuple
    radii: np.ndarray
    z: np.ndarray
    widths: np.ndarray
    membership: np.ndarray
    filament_matrix: np.ndarray
    inductance: np.ndarray
    distribution: np.ndarray

    @property
    def sections(self):
        return self.inductance.shape[0]

    def uniform_inductance(self):
        """Single-turn inductance with uniform current per unit width over all filaments (no equal-flux
        redistribution), for comparison with closed forms."""
        share = self.widths/self.widths.sum()
        return float(share@self.filament_matrix@share)

    def field(self, r, z):
        """Per-section (psi, Br, Bz) per ampere of section current, each shaped like np.broadcast(r, z)."""
        out = []
        for s in range(self.sections):
            psi = br = bz = 0.
            for a, h, d in zip(self.radii, self.z, self.distribution[:, s]):
                if d == 0:
                    continue
                p, b1, b2 = loop_field(r, z, a, h)
                psi, br, bz = psi+d*p, br+d*b1, bz+d*b2
            out.append((psi, br, bz))
        return out

    def sources(self, waveforms, cutoff=1e-9):
        """Coil filaments driven by per-section current waveforms (amperes): filament f carries
        sum_s D[f, s] I_s(t). Returns Coil objects grouped by waveform for VacuumField."""
        coils = []
        for s, wave in enumerate(waveforms):
            scale = np.abs(self.distribution[:, s]).max()
            for k, (a, h, d) in enumerate(zip(self.radii, self.z, self.distribution[:, s])):
                if abs(d) > cutoff*scale:
                    coils.append(Coil(a, h, wave, turns=float(abs(d)), polarity=int(np.sign(d)), name=f"s{s}f{k}"))
        return tuple(coils)


def section_coil(straps, filaments_per_strap=8):
    """Discretize straps into equal-width filaments and reduce to sections.

    Filaments of one (section, turn) link equal flux (perfect conductor on the circuit time scale):
    turn inductance L_t = (G^T M^-1 G)^-1 for the filament-to-turn incidence G. Series turns give the
    section matrix L = H^T L_t H for the turn-to-section incidence H, and filament currents per unit
    section current D = M^-1 G L_t H."""
    straps = tuple(straps)
    turns = sorted({(s.section, s.turn) for s in straps})
    index = {key: k for k, key in enumerate(turns)}
    radii, z, widths, member, turn_of = [], [], [], [], []
    for strap in straps:
        d = strap.width/filaments_per_strap
        for j in range(filaments_per_strap):
            radii.append(strap.radius); z.append(strap.z+(j+.5)*d-strap.width/2); widths.append(d)
            member.append(strap.section); turn_of.append(index[(strap.section, strap.turn)])
    radii, z, widths, member, turn_of = (np.array(v) for v in (radii, z, widths, member, turn_of))
    n_s = int(member.max())+1
    if set(member) != set(range(n_s)):
        raise ValueError("Sections must be numbered 0..n-1 without gaps")
    m = filament_inductance(radii, z, widths)
    g = np.zeros((member.size, len(turns)))
    g[np.arange(member.size), turn_of] = 1.
    h = np.zeros((len(turns), n_s))
    h[np.arange(len(turns)), [sec for sec, _ in turns]] = 1.
    m_inv_g = np.linalg.solve(m, g)
    l_turn = np.linalg.inv(g.T@m_inv_g)
    l_turn = .5*(l_turn+l_turn.T)
    l_sec = h.T@l_turn@h
    return SectionCoil(straps, radii, z, widths, member, m, .5*(l_sec+l_sec.T), m_inv_g@l_turn@h)


@dataclass(frozen=True)
class Crowbar:
    """Switch across the section terminals: armed `delay` after the module capacitor voltage first crosses
    zero (close to current peak), or at absolute `time` if given, and closed then.
    unidirectional (diode-stack crowbar, needs inductance > 0): conducts only the coil's freewheeling current;
    it opens when that current falls to zero (e.g. a later bank drives the coil harder, or a neighbouring
    coil's flux pushes back) and re-closes whenever the section voltage turns forward-biased again."""
    inductance: float = 0.
    resistance: float = 0.
    delay: float = 0.
    time: float | None = None
    unidirectional: bool = False


@dataclass(frozen=True)
class Module:
    """Capacitor C charged to V, in series with stray L (capacitor, switch, feed) and R, switched onto a
    section's terminals at fire_time. Sign of V sets the field direction (bias banks are negative).
    unidirectional: the start switch stops conducting at the first current zero (ignitron, diode-backed
    gap); after a crowbar this isolates the ringing capacitor. Otherwise it conducts both ways."""
    capacitance: float
    voltage: float
    inductance: float
    resistance: float = 0.
    section: int = 0
    fire_time: float = 0.
    crowbar: Crowbar | None = None
    unidirectional: bool = False
    name: str = "module"

    def __post_init__(self):
        if not (self.capacitance > 0 and self.inductance > 0 and self.resistance >= 0 and np.isfinite(self.voltage)):
            raise ValueError("Module needs C > 0, stray L > 0 (parallel ideal sources are singular), R >= 0")


@dataclass
class CircuitResult:
    time: np.ndarray                 # [s]
    module_current: np.ndarray       # (n_t, n_m) [A], positive out of the capacitor into the coil
    capacitor_voltage: np.ndarray    # (n_t, n_m) [V]
    crowbar_current: np.ndarray      # (n_t, n_m) [A], diverted from the section back to the return
    section_current: np.ndarray      # (n_t, n_s) [A]
    section_voltage: np.ndarray      # (n_t, n_s) [V]
    dissipated: np.ndarray           # (n_t,) [J] in all resistances
    plasma_work: np.ndarray          # (n_t,) [J] delivered through the prescribed plasma flux linkage
    crowbar_times: tuple             # closing time per module (None if never)
    modules: tuple
    open_times: tuple = ()           # unidirectional start switches: time each opened (None if never)
    release_times: tuple = ()        # one-way crowbars: first time each opened (None if never)
    crowbar_log: tuple = ()          # per module: crowbar transitions (time, closed)

    def energy_residual(self, inductance):
        """Initial capacitor energy minus (capacitor + magnetic + dissipated + plasma work), per time [J]."""
        c = np.array([m.capacitance for m in self.modules])
        ls = np.array([m.inductance for m in self.modules])
        lc = np.array([m.crowbar.inductance if m.crowbar else 0. for m in self.modules])
        e_cap = .5*(c*self.capacitor_voltage**2).sum(axis=1)
        e_mag = (.5*np.einsum("ti,ij,tj->t", self.section_current, inductance, self.section_current)
                 +.5*(ls*self.module_current**2).sum(axis=1)+.5*(lc*self.crowbar_current**2).sum(axis=1))
        return e_cap[0]-(e_cap+e_mag+self.dissipated+self.plasma_work)

    def section_waveforms(self):
        """Smooth per-section current waveforms (Callback, amperes) for the vacuum-field solver."""
        out = []
        t = self.time
        for s in range(self.section_current.shape[1]):
            spline = CubicSpline(t, self.section_current[:, s])
            rate = spline.derivative()
            out.append(Callback(lambda x, f=spline, t0=t[0]: np.where(np.asarray(x) < t0, f(t0), f(np.maximum(x, t0))),
                                lambda x, f=rate, t0=t[0]: np.where(np.asarray(x) < t0, 0., f(np.maximum(x, t0))),
                                f"circuit section {s} current [A]"))
        return out


def simulate(inductance, modules, t_end, coil_resistance=0., plasma_flux=None, samples=4000, rtol=1e-10):
    """Integrate the switched linear circuit from the earliest fire time to t_end.

    Branch equations (s = module's section):
      module (switch closed):  L_m dI_m/dt + R_m I_m = V_C,m - V_s,     C_m dV_C,m/dt = -I_m
      crowbar (closed):        L_c dI_c/dt + R_c I_c = V_s
      section:                 V_s = sum_s' L_ss' dI_s'/dt + R_coil,s I_s + dPhi_p,s/dt,
                               I_s = sum(I_m) - sum(I_c) over the section's branches.
    plasma_flux(t) -> (Phi_p, dPhi_p/dt), arrays per section [Wb], is an optional prescribed linkage.
    samples: number of uniform output times, or an array of output times (crowbar closings are added)."""
    modules = tuple(modules)
    l_sec = np.atleast_2d(np.asarray(inductance, dtype=float))
    n_s, n_m = l_sec.shape[0], len(modules)
    r_coil = np.broadcast_to(np.asarray(coil_resistance, dtype=float), (n_s,))
    sec = np.array([m.section for m in modules])
    if sec.max() >= n_s:
        raise ValueError("Module attached to a nonexistent section")
    lm = np.array([m.inductance for m in modules]); rm = np.array([m.resistance for m in modules])
    cm = np.array([m.capacitance for m in modules])
    lc = np.array([m.crowbar.inductance if m.crowbar else 0. for m in modules])
    rc = np.array([m.crowbar.resistance if m.crowbar else 0. for m in modules])
    if any(m.crowbar and m.crowbar.unidirectional and m.crowbar.inductance <= 0 for m in modules):
        raise ValueError("A one-way crowbar needs inductance > 0")
    fired = np.zeros(n_m, bool); shorted = np.zeros(n_m, bool); opened = np.zeros(n_m, bool)
    armed = np.zeros(n_m, bool)
    open_at = [None]*n_m
    closing = [m.crowbar.time if (m.crowbar and m.crowbar.time is not None) else None for m in modules]
    log = [[] for _ in range(n_m)]        # crowbar transitions (time, closed)
    plasma = plasma_flux or (lambda t: (np.zeros(n_s), np.zeros(n_s)))
    sign = np.array([np.sign(m.voltage) or 1. for m in modules])

    def section_currents(y):
        i_m, i_c = y[:n_m]*(fired & ~opened), y[2*n_m:3*n_m]*shorted
        return np.bincount(sec, weights=i_m-i_c, minlength=n_s)

    def rhs(t, y):
        i_m, v_c, i_c = y[:n_m], y[n_m:2*n_m], y[2*n_m:3*n_m]
        am, ac = np.nonzero(fired & ~opened)[0], np.nonzero(shorted)[0]
        na, nc = am.size, ac.size
        size = na+nc+n_s
        a = np.zeros((size, size)); b = np.zeros(size)
        i_s = section_currents(y)
        _, dphi = plasma(t)
        for row, m in enumerate(am):          # L dI_m + V_s = V_C - R I_m
            a[row, row] = lm[m]; a[row, na+nc+sec[m]] = 1.; b[row] = v_c[m]-rm[m]*i_m[m]
        for k, m in enumerate(ac):            # L_c dI_c - V_s = -R_c I_c
            row = na+k
            a[row, row] = lc[m]; a[row, na+nc+sec[m]] = -1.; b[row] = -rc[m]*i_c[m]
        for s in range(n_s):                  # V_s - sum L_ss' dI_s' = R I_s + dPhi_p
            row = na+nc+s
            a[row, row] = 1.
            for col, m in enumerate(am):
                a[row, col] -= l_sec[s, sec[m]]
            for k, m in enumerate(ac):
                a[row, na+k] += l_sec[s, sec[m]]
            b[row] = r_coil[s]*i_s[s]+dphi[s]
        x = np.linalg.solve(a, b) if size else np.zeros(0)
        dy = np.zeros_like(y)
        dy[am] = x[:na]
        dy[n_m+am] = -i_m[am]/cm[am]
        dy[2*n_m+ac] = x[na:na+nc]
        v_s = x[na+nc:]
        power = (rm*i_m**2*(fired & ~opened)).sum()+(rc*i_c**2*shorted).sum()+(r_coil*i_s**2).sum()
        dy[3*n_m] = power
        dy[3*n_m+1] = (i_s*dphi).sum()
        return dy, v_s

    # Switch quantities that act at a zero crossing (from + to -): crowbar trigger (capacitor voltage), start
    # switch opening (module current), one-way crowbar release (freewheeling current) and re-close (voltage).
    def quantity(kind, k, tt, yy):
        if kind == "trigger":
            return yy[n_m+k]*sign[k]
        if kind == "gate":
            return yy[k]*sign[k]
        if kind == "release":
            return -yy[2*n_m+k]*sign[k]
        return rhs(tt, yy)[1][sec[k]]*sign[k]   # reclose: section voltage forward-biases the crowbar

    def watched():
        out = []
        for k, m in enumerate(modules):
            if m.crowbar and m.crowbar.time is None and fired[k] and closing[k] is None:
                out.append(("trigger", k))
            if m.unidirectional and fired[k] and not opened[k]:
                out.append(("gate", k))
            if m.crowbar and m.crowbar.unidirectional and shorted[k]:
                out.append(("release", k))
            if m.crowbar and m.crowbar.unidirectional and armed[k] and not shorted[k]:
                out.append(("reclose", k))
        return out

    def tolerance(kind, k, yy):
        if kind == "trigger" or kind == "reclose":
            return 1e-9*abs(modules[k].voltage)
        return 1e-9*max(1., np.abs(yy[:3*n_m]).max())

    def apply(kind, k, tt, yy):
        if kind == "trigger":
            closing[k] = tt+modules[k].crowbar.delay
        elif kind == "gate":
            opened[k] = True; open_at[k] = tt; yy[k] = 0.
        elif kind == "release":
            shorted[k] = False; yy[2*n_m+k] = 0.; log[k].append((tt, False))
        else:
            shorted[k] = True; yy[2*n_m+k] = 0.; log[k].append((tt, True))

    y = np.zeros(3*n_m+2)
    y[n_m:2*n_m] = [m.voltage for m in modules]
    t = min(m.fire_time for m in modules)

    def schedule():
        pending = [m.fire_time for k, m in enumerate(modules) if not fired[k] and m.fire_time > t]
        pending += [c for k, c in enumerate(closing) if c is not None and not armed[k] and c > t]
        return min(pending+[t_end])

    def switch(now, yy):
        for k, m in enumerate(modules):
            if not fired[k] and m.fire_time <= now:
                fired[k] = True
            if closing[k] is not None and not armed[k] and closing[k] <= now and fired[k]:
                armed[k] = True; shorted[k] = True; yy[2*n_m+k] = 0.; log[k].append((now, True))

    switch(t, y)
    t_start, segments = t, []
    while t < t_end:
        stop = schedule()
        kinds = watched()
        events = []
        for kind, k in kinds:
            def event(tt, yy, kind=kind, k=k):
                return quantity(kind, k, tt, yy)
            event.terminal, event.direction = True, -1
            events.append(event)
        scale = max(1., np.abs(y[:3*n_m]).max(), max(abs(m.voltage) for m in modules))
        sol = solve_ivp(lambda tt, yy: rhs(tt, yy)[0], (t, stop), y, method="DOP853", rtol=rtol,
                        atol=1e-13*scale, events=events or None, dense_output=True)
        if sol.status < 0:
            raise RuntimeError(sol.message)
        segments.append((sol.t[0], sol.t[-1], sol.sol))
        t, y = sol.t[-1], sol.y[:, -1].copy()
        if sol.status == 1:
            # Apply every switch quantity that has reached zero at the stop time, not only the one that stopped
            # the integration: symmetric modules cross together and would otherwise be missed for good.
            due = [(kind, k) for (kind, k), hits in zip(kinds, sol.t_events)
                   if hits.size or quantity(kind, k, t, y) <= tolerance(kind, k, y)]
            for kind, k in due:
                apply(kind, k, t, y)
        switch(t, y)
    if np.ndim(samples):
        grid = np.asarray(samples, dtype=float)
        grid = grid[(grid >= t_start) & (grid <= t_end)]
    else:
        grid = np.linspace(t_start, t_end, samples)
    grid = np.union1d(grid, [tt for entries in log for tt, _ in entries if t_start <= tt <= t_end])
    y_grid = np.empty((grid.size, y.size))
    for n, tt in enumerate(grid):
        a, b, dense = next(s for s in segments if s[0] <= tt <= s[1]) if tt < segments[-1][1] else segments[-1]
        y_grid[n] = dense(min(max(tt, a), b))
    # Replay switch states on the sample grid for derived quantities.
    i_m = np.zeros((grid.size, n_m)); i_c = np.zeros((grid.size, n_m)); i_s = np.zeros((grid.size, n_s))
    v_s = np.zeros((grid.size, n_s))
    for n, tt in enumerate(grid):
        fired[:] = [m.fire_time <= tt for m in modules]
        shorted[:] = [bool([closed for when, closed in log[k] if when <= tt][-1:] == [True]) for k in range(n_m)]
        opened[:] = [open_at[k] is not None and open_at[k] <= tt for k in range(n_m)]
        i_m[n] = y_grid[n, :n_m]*(fired & ~opened); i_c[n] = y_grid[n, 2*n_m:3*n_m]*shorted
        i_s[n] = section_currents(y_grid[n])
        v_s[n] = rhs(tt, y_grid[n])[1]
    first_close = tuple(next((when for when, closed in log[k] if closed), None) for k in range(n_m))
    releases = tuple(next((when for when, closed in log[k] if not closed), None) for k in range(n_m))
    return CircuitResult(grid, i_m, y_grid[:, n_m:2*n_m], i_c, i_s, v_s, y_grid[:, 3*n_m], y_grid[:, 3*n_m+1],
                         first_close, modules, tuple(open_at), releases, tuple(tuple(e) for e in log))


class LinearWaveform:
    """Mutable current waveform I(t) = current + rate (t - start), updated by a coupled circuit each step."""
    def __init__(self, current=0., rate=0., start=0.):
        self.current, self.rate, self.start = current, rate, start

    def value(self, t):
        return self.current+self.rate*(np.asarray(t, dtype=float)-self.start)

    def derivative(self, t):
        return np.zeros_like(np.asarray(t, dtype=float))+self.rate


class CircuitStepper:
    """The circuit of `simulate`, advanced in lock-step with a plasma solver whose coil is a flux conserver.

    The solver's coil boundary carries an image current I_img,s per section on top of the prescribed
    vacuum current I_vac,s, and the coil links the vacuum flux L I_vac only. KCL at each section node:
    I_vac,s = sum(I_m) - sum(I_c) - I_img,s, so the image current is a load the banks must supply and the
    stray inductances limit. Each solver step: RK4 over [t, t+dt] with I_img linear in time; switch events
    (fire, crowbar arming at the capacitor-voltage zero, start-switch opening, one-way crowbar release and
    re-close) are located inside the step, as in `simulate`."""
    def __init__(self, inductance, modules, coil_resistance=0., start=None):
        self.l = np.atleast_2d(np.asarray(inductance, dtype=float))
        self.modules = tuple(modules)
        n_s, n_m = self.l.shape[0], len(self.modules)
        self.n_s, self.n_m = n_s, n_m
        if any(m.crowbar and m.crowbar.unidirectional and m.crowbar.inductance <= 0 for m in self.modules):
            raise ValueError("A one-way crowbar needs inductance > 0")
        self.r_coil = np.broadcast_to(np.asarray(coil_resistance, dtype=float), (n_s,))
        self.sec = np.array([m.section for m in self.modules])
        self.lm = np.array([m.inductance for m in self.modules]); self.rm = np.array([m.resistance for m in self.modules])
        self.cm = np.array([m.capacitance for m in self.modules])
        self.lc = np.array([m.crowbar.inductance if m.crowbar else 0. for m in self.modules])
        self.rc = np.array([m.crowbar.resistance if m.crowbar else 0. for m in self.modules])
        self.sign = np.array([np.sign(m.voltage) or 1. for m in self.modules])
        self.t = min(m.fire_time for m in self.modules) if start is None else start
        self.y = np.zeros(3*n_m+2)
        self.y[n_m:2*n_m] = [m.voltage for m in self.modules]
        self.fired = np.array([m.fire_time <= self.t for m in self.modules])
        self.shorted = np.zeros(n_m, bool); self.opened = np.zeros(n_m, bool); self.armed = np.zeros(n_m, bool)
        self.closing = [m.crowbar.time if (m.crowbar and m.crowbar.time is not None) else None for m in self.modules]
        self.open_at = [None]*n_m
        self.log = [[] for _ in range(n_m)]   # crowbar transitions (time, closed)
        self.image = np.zeros(n_s)
        self.waveforms = [LinearWaveform(start=self.t) for _ in range(n_s)]

    def _solve(self, t, y, image_rate):
        n_m, n_s = self.n_m, self.n_s
        active = self.fired & ~self.opened
        am, ac = np.nonzero(active)[0], np.nonzero(self.shorted)[0]
        na, nc = am.size, ac.size
        a = np.zeros((na+nc+n_s, na+nc+n_s)); b = np.zeros(na+nc+n_s)
        i_m, v_c, i_c = y[:n_m], y[n_m:2*n_m], y[2*n_m:3*n_m]
        i_vac = self.vacuum_current(y)
        for row, m in enumerate(am):
            a[row, row] = self.lm[m]; a[row, na+nc+self.sec[m]] = 1.; b[row] = v_c[m]-self.rm[m]*i_m[m]
        for k, m in enumerate(ac):
            row = na+k
            a[row, row] = self.lc[m]; a[row, na+nc+self.sec[m]] = -1.; b[row] = -self.rc[m]*i_c[m]
        for s in range(n_s):  # V_s - L (sum dI_m - sum dI_c) = R I_vac - L dI_img
            row = na+nc+s
            a[row, row] = 1.
            for col, m in enumerate(am):
                a[row, col] -= self.l[s, self.sec[m]]
            for k, m in enumerate(ac):
                a[row, na+k] += self.l[s, self.sec[m]]
            b[row] = self.r_coil[s]*i_vac[s]-self.l[s]@image_rate
        x = np.linalg.solve(a, b) if a.size else np.zeros(0)
        return am, ac, x, i_vac

    def _rhs(self, t, y, image_rate):
        n_m, n_s = self.n_m, self.n_s
        am, ac, x, i_vac = self._solve(t, y, image_rate)
        na, nc = am.size, ac.size
        active = self.fired & ~self.opened
        i_m, i_c = y[:n_m], y[2*n_m:3*n_m]
        dy = np.zeros_like(y)
        dy[am] = x[:na]; dy[n_m+am] = -i_m[am]/self.cm[am]; dy[2*n_m+ac] = x[na:na+nc]
        dy[3*n_m] = (self.rm*i_m**2*active).sum()+(self.rc*i_c**2*self.shorted).sum()+(self.r_coil*i_vac**2).sum()
        rate_vac = np.bincount(self.sec, weights=dy[:n_m]*active-dy[2*n_m:3*n_m]*self.shorted, minlength=n_s)-image_rate
        dy[3*n_m+1] = float(self.image@(self.l@rate_vac))   # V_s I_img: power into the plasma through the image load
        return dy

    def vacuum_current(self, y=None, image=None):
        y = self.y if y is None else y
        image = self.image if image is None else image
        active = self.fired & ~self.opened
        branch = np.bincount(self.sec, weights=y[:self.n_m]*active-y[2*self.n_m:3*self.n_m]*self.shorted,
                             minlength=self.n_s)
        return branch-image

    def branch(self, image_rate=None):
        """Section branch currents sum(I_m) - sum(I_c) = I_vac + I_img [A] at self.t, and their rates [A/s] with the
        image changing at `image_rate` (the current a directly driven load, e.g. a rod's boundary, carries)."""
        image_rate = np.zeros(self.n_s) if image_rate is None else image_rate
        dy = self._rhs(self.t, self.y, image_rate)
        active = self.fired & ~self.opened
        n_m = self.n_m
        current = np.bincount(self.sec, weights=self.y[:n_m]*active-self.y[2*n_m:3*n_m]*self.shorted, minlength=self.n_s)
        rate = np.bincount(self.sec, weights=dy[:n_m]*active-dy[2*n_m:3*n_m]*self.shorted, minlength=self.n_s)
        return current, rate

    def connected(self):
        """Sections with a closed path for net current: a conducting module or a closed crowbar."""
        active = self.fired & ~self.opened
        out = np.zeros(self.n_s, bool)
        out[self.sec[active | self.shorted]] = True
        return out

    def vacuum_rate(self, image_rate):
        return self._rhs(self.t, self.y, image_rate)[:self.n_m]

    def _watched(self):
        out = []
        for k, m in enumerate(self.modules):
            if m.crowbar and m.crowbar.time is None and self.fired[k] and self.closing[k] is None:
                out.append(("trigger", k))
            if m.unidirectional and self.fired[k] and not self.opened[k]:
                out.append(("gate", k))
            if m.crowbar and m.crowbar.unidirectional and self.shorted[k]:
                out.append(("release", k))
            if m.crowbar and m.crowbar.unidirectional and self.armed[k] and not self.shorted[k]:
                out.append(("reclose", k))
        return out

    def _quantity(self, kind, k, t, y, image, image_rate):
        """Switch quantity that acts on a + to - zero crossing (see simulate)."""
        n_m = self.n_m
        if kind == "trigger":
            return y[n_m+k]*self.sign[k]
        if kind == "gate":
            return y[k]*self.sign[k]
        if kind == "release":
            return -y[2*n_m+k]*self.sign[k]
        saved = self.image
        self.image = image
        am, ac, x, _ = self._solve(t, y, image_rate)
        self.image = saved
        return x[am.size+ac.size+self.sec[k]]*self.sign[k]

    def _tolerance(self, kind, k, y):
        if kind in ("trigger", "reclose"):
            return 1e-9*abs(self.modules[k].voltage)
        return 1e-9*max(1., np.abs(y[:3*self.n_m]).max())

    def _apply(self, kind, k):
        n_m = self.n_m
        if kind == "trigger":
            self.closing[k] = self.t+self.modules[k].crowbar.delay
        elif kind == "gate":
            self.opened[k] = True; self.open_at[k] = self.t; self.y[k] = 0.
        elif kind == "release":
            self.shorted[k] = False; self.y[2*n_m+k] = 0.; self.log[k].append((self.t, False))
        else:
            self.shorted[k] = True; self.y[2*n_m+k] = 0.; self.log[k].append((self.t, True))

    def advance(self, t_new, image_new):
        """Advance from self.t to t_new with the image current linear between its old and new values.
        Sub-steps (RK4) end at fire and crowbar-arming times and at switch events located by linear
        interpolation of their quantities; every quantity at zero then is applied together. The per-section
        LinearWaveforms are then refreshed (value and slope at t_new)."""
        if t_new <= self.t:
            return
        image_new = np.asarray(image_new, dtype=float)
        image_rate = (image_new-self.image)/(t_new-self.t)
        n_m = self.n_m
        start = self.t
        while self.t < t_new:
            fires = [m.fire_time for k, m in enumerate(self.modules) if not self.fired[k] and self.t < m.fire_time < t_new]
            arms = [c for k, c in enumerate(self.closing) if c is not None and not self.armed[k] and self.t < c < t_new]
            t_end = min(fires+arms+[t_new])
            image_end = self.image+image_rate*(t_end-self.t)
            y1 = self._rk4(t_end-self.t, image_rate)
            watched = self._watched()
            first = 1.
            for kind, k in watched:
                g0 = self._quantity(kind, k, self.t, self.y, self.image, image_rate)
                g1 = self._quantity(kind, k, t_end, y1, image_end, image_rate)
                if g0 > 0 >= g1 and g0/(g0-g1) < first:
                    first = g0/(g0-g1)
            if first < 1.:
                t_end = self.t+first*(t_end-self.t)
                image_end = self.image+image_rate*(t_end-self.t)
                y1 = self._rk4(t_end-self.t, image_rate)
            self.y, self.t, self.image = y1, t_end, image_end
            if first < 1.:
                due = [(kind, k) for kind, k in watched
                       if self._quantity(kind, k, self.t, self.y, self.image, image_rate) <= self._tolerance(kind, k, self.y)]
                for kind, k in due:
                    self._apply(kind, k)
            for k, m in enumerate(self.modules):
                if not self.fired[k] and m.fire_time <= self.t:
                    self.fired[k] = True
                if self.closing[k] is not None and not self.armed[k] and self.closing[k] <= self.t and self.fired[k]:
                    self.armed[k] = True; self.shorted[k] = True; self.y[2*n_m+k] = 0.; self.log[k].append((self.t, True))
        self.image = image_new
        self.refresh(image_rate, relax=t_new-start)

    def _rk4(self, dt, image_rate):
        y0, image0, t0 = self.y, self.image, self.t

        def f(tt, yy):
            self.image = image0+image_rate*(tt-t0)
            return self._rhs(tt, yy, image_rate)
        k1 = f(t0, y0); k2 = f(t0+dt/2, y0+dt/2*k1); k3 = f(t0+dt/2, y0+dt/2*k2); k4 = f(t0+dt, y0+dt*k3)
        self.image = image0
        return y0+dt/6*(k1+2*k2+2*k3+k4)

    def refresh(self, image_rate=None, relax=None):
        """Set each section's source waveform for the next solver step. Without `relax` the waveform takes the
        circuit's vacuum current and slope at self.t. With `relax` (the last step length) it stays continuous:
        it keeps its own value at self.t and catches up with the circuit over one further step,
        slope = circuit slope + (circuit current - waveform value)/relax. A source that jumped at a step
        boundary would change the field with no matching work (an open section's current is set algebraically
        by its image load, so it can appear suddenly)."""
        image_rate = np.zeros(self.n_s) if image_rate is None else image_rate
        current = self.vacuum_current()
        dy = self._rhs(self.t, self.y, image_rate)
        active = self.fired & ~self.opened
        rate = np.bincount(self.sec, weights=dy[:self.n_m]*active-dy[2*self.n_m:3*self.n_m]*self.shorted,
                           minlength=self.n_s)-image_rate
        for s, wave in enumerate(self.waveforms):
            if relax is None or relax <= 0:
                wave.current, wave.rate, wave.start = float(current[s]), float(rate[s]), self.t
            else:
                held = float(wave.value(self.t))
                wave.current, wave.rate, wave.start = held, float(rate[s]+(current[s]-held)/relax), self.t

    def run_to(self, t_end, dt):
        """Advance without plasma (image held) in steps of at most dt, e.g. through the bias phase; the
        waveforms then take the circuit state exactly (no plasma has seen them yet)."""
        while self.t < t_end-1e-15:
            self.advance(min(self.t+dt, t_end), self.image)
        self.refresh()

    def energy_residual(self):
        """Initial capacitor energy minus capacitor, vacuum-coil, stray and crowbar magnetic energy,
        dissipation and the work delivered through the image load [J]."""
        n_m = self.n_m
        i_vac = self.vacuum_current()
        i_m = self.y[:n_m]*(self.fired & ~self.opened); i_c = self.y[2*n_m:3*n_m]*self.shorted
        e0 = sum(.5*m.capacitance*m.voltage**2 for m in self.modules)
        return e0-(.5*(self.cm*self.y[n_m:2*n_m]**2).sum()+.5*i_vac@self.l@i_vac+.5*(self.lm*i_m**2).sum()
                   +.5*(self.lc*i_c**2).sum()+self.y[3*n_m]+self.y[3*n_m+1])

    def snapshot(self):
        n_m = self.n_m
        return dict(time_s=self.t, vacuum_current_a=self.vacuum_current().tolist(), image_current_a=self.image.tolist(),
                    module_current_a=(self.y[:n_m]*(self.fired & ~self.opened)).tolist(),
                    capacitor_voltage_v=self.y[n_m:2*n_m].tolist(), crowbar_current_a=(self.y[2*n_m:3*n_m]*self.shorted).tolist(),
                    dissipated_j=float(self.y[3*n_m]), image_work_j=float(self.y[3*n_m+1]))
