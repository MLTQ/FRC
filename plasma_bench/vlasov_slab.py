"""Conservative 1D1V collisionless electrostatic slab with physical species masses."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import epsilon_0, elementary_charge as EV
from scipy.special import ndtr


def gaussian_interval(lo, hi):
    """Standard Gaussian probability; avoid cancellation in the positive tail."""
    lo, hi = np.broadcast_arrays(lo, hi)
    return np.where(lo >= 0, ndtr(-lo)-ndtr(-hi), ndtr(hi)-ndtr(lo))


def advect(f, courant, axis):
    """Constant-along-line transport, zero inflow, minmod linear reconstruction."""
    g = np.moveaxis(f, axis, -1)
    c = np.broadcast_to(courant, g.shape[:-1])[..., None]
    if not np.isfinite(c).all() or np.max(abs(c)) > .8:
        raise ValueError('Requested advection exceeds the explicit .8 CFL bound')
    if not np.isfinite(g).all() or np.min(g) < 0:
        raise ValueError('Nonfinite or negative phase density; no repair')
    left = np.empty_like(g); right = np.empty_like(g)
    left[..., 0] = g[..., 0]; left[..., 1:] = np.diff(g, axis=-1)
    right[..., -1] = -g[..., -1]; right[..., :-1] = left[..., 1:]
    slope = np.where(np.signbit(left) == np.signbit(right),
                     np.sign(left)*np.minimum(abs(left), abs(right)), 0.)
    magnitude, direction = abs(c), np.sign(c)
    outgoing = magnitude*(g+.5*direction*(1-magnitude)*slope)
    result = (1-magnitude)*(g-.5*magnitude*direction*slope)
    result[..., 1:] += np.where(c >= 0, outgoing[..., :-1], 0.)
    result[..., :-1] += np.where(c < 0, outgoing[..., 1:], 0.)
    losses = np.stack((np.where(c[..., 0] < 0, outgoing[..., 0], 0.),
                       np.where(c[..., 0] > 0, outgoing[..., -1], 0.)), axis=-1)
    if not np.isfinite(result).all() or np.min(result) < 0:
        raise ValueError('Invalid advected phase density; no clipping')
    return np.moveaxis(result, -1, axis), losses


@dataclass(frozen=True)
class KineticSpecies:
    name: str
    mass: float
    charge: float
    temperature_ev: float
    velocity_cells: int
    velocity_sigma: float = 8.

    @property
    def thermal_speed(self):
        return np.sqrt(EV*self.temperature_ev/self.mass)

    @property
    def s_edges(self):
        return np.linspace(-self.velocity_sigma, self.velocity_sigma, self.velocity_cells+1)

    @property
    def s(self):
        return (self.s_edges[1:]+self.s_edges[:-1])/2

    @property
    def ds(self):
        return 2*self.velocity_sigma/self.velocity_cells

    def raw(self, f):
        """Exact velocity moments of the represented piecewise-constant density."""
        s, d = self.s, self.ds
        weights = np.stack((np.ones_like(s), s, s*s+d*d/12, s**3+s*d*d/4))
        return (f @ weights.T)*d

    def profiles(self, f, density):
        raw = self.raw(f); n, first, second, third = raw.T
        u = np.divide(first, n, out=np.zeros_like(first), where=n > 0)
        variance = second-2*u*first+u*u*n
        central_third = third-3*u*second+3*u*u*first-u**3*n
        return dict(number_m3=density*n, velocity_m_s=self.thermal_speed*u,
                    current_a_m2=self.charge*density*self.thermal_speed*first,
                    parallel_pressure_pa=EV*self.temperature_ev*density*variance,
                    perpendicular_pressure_pa=EV*self.temperature_ev*density*n,
                    parallel_heat_w_m2=.5*self.mass*density*self.thermal_speed**3*central_third,
                    raw=raw)


class VlasovSlab:
    def __init__(self, species, spatial_cells, self_consistent=True,
                 density=1e14, half_width=.01, ends=(-.2, .2)):
        if spatial_cells < 2 or density <= 0 or not ends[0] < -half_width < half_width < ends[1]:
            raise ValueError('Invalid slab geometry or density')
        if any(s.mass <= 0 or s.temperature_ev <= 0 or s.velocity_cells < 2 or s.velocity_sigma <= 0 for s in species):
            raise ValueError('Invalid kinetic species')
        self.species = tuple(species); self.edges = np.linspace(*ends, spatial_cells+1)
        self.dx = (ends[1]-ends[0])/spatial_cells
        self.density, self.half_width = density, half_width
        self.self_consistent = self_consistent

    def initial(self):
        occupied = np.maximum(0., np.minimum(self.edges[1:], self.half_width)
                              -np.maximum(self.edges[:-1], -self.half_width))/self.dx
        return tuple(occupied[:, None]*gaussian_interval(s.s_edges[:-1], s.s_edges[1:])[None, :]/s.ds
                     for s in self.species)

    def electric_field(self, state):
        charge = sum(s.charge*self.density*np.sum(f, axis=1)*s.ds for s, f in zip(self.species, state))
        return np.r_[0., np.cumsum(charge)*self.dx/epsilon_0]

    def timestep(self, state, cfl=.65):
        field = self.electric_field(state) if self.self_consistent else np.zeros(1)
        limits = [cfl*self.dx/(s.velocity_sigma*s.thermal_speed) for s in self.species]
        for s in self.species:
            accel = abs(s.charge)*np.max(abs(field))/s.mass
            if accel > 0:
                limits.append(cfl*s.ds*s.thermal_speed/accel)
        return min(limits)

    def step(self, state, dt):
        """Strang stream/accelerate/stream; returned losses are species x (N,P,E)."""
        half, xloss, vloss = [], [], []
        for s, f in zip(self.species, state):
            g, loss = advect(f, s.s*s.thermal_speed*dt/(2*self.dx), axis=0)
            half.append(g); xloss.append(self._x_loss(s, loss)); vloss.append(np.zeros(3))
        field = self.electric_field(half)
        centre_field = (field[:-1]+field[1:])/2
        final = []
        for i, (s, g) in enumerate(zip(self.species, half)):
            if self.self_consistent:
                c = s.charge*centre_field*dt/(s.mass*s.thermal_speed*s.ds)
                g, loss = advect(g, c, axis=1)
                vloss[i] = self._v_loss(s, loss)
            g, loss = advect(g, s.s*s.thermal_speed*dt/(2*self.dx), axis=0)
            final.append(g); xloss[i] += self._x_loss(s, loss)
        impulse = .5*epsilon_0*(field[-1]**2-field[0]**2)*dt if self.self_consistent else 0.
        return tuple(final), np.asarray(xloss), np.asarray(vloss), impulse

    def _x_loss(self, s, loss):
        number = self.density*self.dx*s.ds*loss.sum(axis=1)
        v = s.s*s.thermal_speed
        energy = .5*s.mass*s.thermal_speed**2*(s.s**2+s.ds**2/12)+EV*s.temperature_ev
        return np.array([number.sum(), number @ (s.mass*v), number @ energy])

    def _v_loss(self, s, loss):
        number = self.density*self.dx*s.ds*loss.sum(axis=0)
        v = s.thermal_speed*s.s_edges[[0, -1]]
        energy = .5*s.mass*v*v+EV*s.temperature_ev
        return np.array([number.sum(), number @ (s.mass*v), number @ energy])

    def budget(self, state):
        totals = []
        for s, f in zip(self.species, state):
            moments = s.raw(f).sum(axis=0)*self.density*self.dx
            totals.append([moments[0], s.mass*s.thermal_speed*moments[1],
                           EV*s.temperature_ev*(.5*moments[2]+moments[0])])
        e = self.electric_field(state)
        field_energy = epsilon_0*self.dx*np.sum(e[:-1]**2+e[:-1]*e[1:]+e[1:]**2)/6
        totals = np.asarray(totals)
        energy = totals[:, 2].sum()+(field_energy if self.self_consistent else 0.)
        return totals, energy, field_energy
