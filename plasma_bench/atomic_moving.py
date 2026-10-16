"""Moving species control volumes with exact exterior vacuum and Gauss field work."""
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from .atomic_transport import AtomicTransport
from .species_ale_flux import moving_flux
from .species_reconstruction import reconstruct


class MovingAtomicElectrostatic(AtomicTransport):
    """Axial ALE species amounts, freely tracked fronts, fixed remote field ends."""
    def __init__(self, source, faces, reference_density, outer_bounds, *, left_electric=0., cfl=.2, maximum_plasma_phase=.04, mesh_motion='affine'):
        faces = np.asarray(faces, float)
        super().__init__(source, len(faces)-1, faces[-1]-faces[0], reference_density, cfl=cfl)
        self.bounds = np.asarray(outer_bounds, float)
        if self.bounds.shape != (2,) or not np.isfinite(self.bounds).all() or self.bounds[1] <= self.bounds[0]:
            raise ValueError('Ordered finite remote electric-field boundaries required')
        if not np.isfinite([left_electric, maximum_plasma_phase]).all() or maximum_plasma_phase <= 0:
            raise ValueError('Finite boundary E and positive plasma-phase step required')
        self.left_electric = float(left_electric); self.maximum_plasma_phase = float(maximum_plasma_phase)
        if mesh_motion not in ('affine', 'barycentric'):
            raise ValueError('Choose affine-front or barycentric mesh motion')
        self.mesh_motion = mesh_motion
        self.reconstruction_count = 0; self.limited_reconstruction_count = 0
        self.charge = np.r_[source.network.charge, -1.]
        self.geometry(faces)

    def geometry(self, faces):
        if faces.shape != (self.cells+1,) or not np.isfinite(faces).all():
            raise ValueError('Complete finite moving face coordinates required')
        width = np.diff(faces)
        if np.any(width <= 0) or faces[0] <= self.bounds[0] or faces[-1] >= self.bounds[1]:
            raise ValueError('Moving cells crossed or reached a remote field boundary; no remapping')
        return width, .5*(faces[:-1]+faces[1:])

    def primitive(self, faces, amounts):
        width, _ = self.geometry(faces)
        return super().primitive(amounts/width[:, None])

    def traces(self, faces, amounts):
        width, _ = self.geometry(faces); p = self.primitive(faces, amounts); s = self.source
        density = amounts/width[:, None]
        conserved = np.concatenate((density[:, self.number, None], density[:, s.momentum], density[:, s.energy, None]), axis=2)
        left, right, factor = reconstruct(conserved, faces, self.mass_ratio, self.gamma)
        occupied = p[..., 0] > 0
        self.reconstruction_count += int(occupied.sum())
        self.limited_reconstruction_count += int(np.count_nonzero(occupied&(factor < 1.)))
        return p, left, right

    def mesh_velocity(self, faces, amounts, traces=None):
        width, _ = self.geometry(faces)
        p, left, right = self.traces(faces, amounts) if traces is None else traces
        mass = p[..., 0]*self.mass_ratio
        total = mass.sum(axis=1)
        if np.any(total <= 0):
            raise ValueError('Material control volumes must contain matter; exterior has no cells')
        velocity = self.source.speed_scale*(np.sum(mass*p[..., 4], axis=1)/total+self.frame[2])
        mesh = np.zeros(self.cells+1)
        mesh[1:-1] = (width[1:]*velocity[:-1]+width[:-1]*velocity[1:])/(width[:-1]+width[1:])
        # A reconstructed vacuum trace may correctly have n=p=0, where p/n
        # does not define a sound speed. Pose the boundary Riemann problem
        # between the adjacent occupied cell state and exact exterior vacuum.
        for index, trace, sign in ((0, p[0], -1), (-1, p[-1], 1)):
            occupied = trace[:, 0] > 0
            if not np.any(occupied):
                raise ValueError('No material in the front control volume')
            sound = self.source.speed_scale*np.sqrt(self.gamma*trace[occupied, 1]/trace[occupied, 0]/self.mass_ratio[occupied])
            speeds = self.source.speed_scale*(trace[occupied, 4]+self.frame[2])+sign*2*sound/(self.gamma-1)
            mesh[index] = speeds.min() if sign < 0 else speeds.max()
        if self.mesh_motion == 'affine':
            # Coordinate motion only: the physical species velocities are still
            # independent unknowns, transported through these moving faces.
            fraction = (faces[1:-1]-faces[0])/(faces[-1]-faces[0])
            mesh[1:-1] = mesh[0]+fraction*(mesh[-1]-mesh[0])
        return mesh

    def electric(self, amounts):
        return np.r_[self.left_electric, self.left_electric+EV*self.density/EPS0*np.cumsum(amounts[:, self.number]@self.charge)]

    def timestep(self, faces, amounts):
        width, _ = self.geometry(faces); traces = self.traces(faces, amounts); p = traces[0]
        mesh = self.mesh_velocity(faces, amounts, traces)
        sound = self.source.speed_scale*np.sqrt(self.gamma*np.divide(p[..., 1], p[..., 0], out=np.zeros_like(p[..., 0]), where=p[..., 0] > 0)/self.mass_ratio)
        velocity = self.source.speed_scale*(p[..., 4]+self.frame[2])
        speed = np.maximum(abs(velocity-mesh[:-1, None]), abs(velocity-mesh[1:, None]))+sound
        speed = np.max(np.where(p[..., 0] > 0, speed, 0.), axis=1)
        dt = np.min(np.divide(self.cfl*width, speed, out=np.full_like(width, np.inf), where=speed > 0))
        contraction = -np.diff(mesh)
        dt = min(dt, np.min(np.divide(.3*width, contraction, out=np.full_like(width, np.inf), where=contraction > 0)))
        plasma = np.sqrt(EV**2*self.density/EPS0*np.max(p[..., 0]@(self.charge**2/self.source.mass)))
        return min(dt, self.maximum_plasma_phase/plasma) if plasma else dt

    def derivative(self, faces, amounts):
        width, _ = self.geometry(faces); traces = self.traces(faces, amounts)
        mesh = self.mesh_velocity(faces, amounts, traces); _, left, right = traces; s = self.source
        flux = np.zeros((self.cells+1, len(self.number), 5))
        flux[1:-1] = moving_flux(right[:-1], left[1:], self.mass_ratio, self.frame[2], mesh[1:-1]/s.speed_scale)*s.speed_scale
        # The exact vacuum fan front has zero ALE particle, momentum and energy
        # flux. Front speed is the outward envelope of all occupied species.
        change = -np.diff(flux, axis=0); rate = np.zeros_like(amounts)
        rate[:, self.number] = change[..., 0]; rate[:, s.momentum] = change[..., 1:4]; rate[:, s.energy] = change[..., 4]
        field = self.electric(amounts); a, b = field[:-1], field[1:]
        force = (EV/s.momentum_scale)*self.charge*amounts[:, self.number]*((a+b)/2)[:, None]
        rate[:, s.momentum[:, 2]] += force
        # Integrate E times the actual continuity current in each moving cell.
        weighted_left = (2*a+b)[:, None]; weighted_right = (a+2*b)[:, None]
        work = self.charge/6*(width[:, None]*(weighted_left*flux[:-1, :, 0]+weighted_right*flux[1:, :, 0])
            + amounts[:, self.number]*(weighted_left*mesh[:-1, None]+weighted_right*mesh[1:, None]))
        rate[:, s.energy] += work - (s.momentum_scale/EV)*s.speed_scale*self.frame[2]*force
        field_rate = np.r_[0., EV*self.density/EPS0*np.cumsum(rate[:, self.number]@self.charge)]
        scale = EPS0/(6*EV*self.density)
        power = scale*np.sum(np.diff(mesh)*(a*a+a*b+b*b)+width*((2*a+b)*field_rate[:-1]+(a+2*b)*field_rate[1:]))
        power += 3*scale*(field[0]**2*mesh[0]-field[-1]**2*mesh[-1]
            +2*field[-1]*field_rate[-1]*(self.bounds[1]-faces[-1]))
        self.last_diagnostic = dict(electric_work_residual_ev_m_s=float(work.sum()+power),
            left_front_speed_m_s=float(mesh[0]), right_front_speed_m_s=float(mesh[-1]))
        return mesh, rate

    def step(self, faces, amounts, time, dt, *, reactions=True, heat=None):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite moving step required')

        def material(x, q, start, duration):
            width, centre = self.geometry(x)
            initial = np.zeros(self.cells) if heat is None else np.broadcast_to(heat(centre, start), (self.cells,))
            middle = np.zeros(self.cells) if heat is None else np.broadcast_to(heat(centre, start+duration/2), (self.cells,))
            return width[:, None]*np.asarray([self.source.step(row, duration, self.density,
                initial_heat=initial[j], midpoint_heat=middle[j]) for j, row in enumerate(q/width[:, None])])

        q = material(faces, amounts, time, dt/2) if reactions else amounts
        values = (faces, q); k1 = self.derivative(*values)
        k2 = self.derivative(*(u+.5*dt*k for u, k in zip(values, k1)))
        k3 = self.derivative(*(u+.5*dt*k for u, k in zip(values, k2)))
        k4 = self.derivative(*(u+dt*k for u, k in zip(values, k3)))
        result = tuple(u+dt/6*(a+2*b+2*d+e) for u, a, b, d, e in zip(values, k1, k2, k3, k4))
        self.primitive(*result)
        return (result[0], material(*result, time+dt/2, dt/2)) if reactions else result

    def sample(self, faces, amounts, points):
        """Actual piecewise-constant density state, exactly zero outside support."""
        points = np.asarray(points); indices = np.searchsorted(faces, points, side='right')-1
        result = np.zeros(points.shape+(self.source.size,)); inside = (indices >= 0)&(indices < self.cells)
        result[inside] = (amounts/np.diff(faces)[:, None])[indices[inside]]
        return result

    def totals(self, faces, amounts):
        width, _ = self.geometry(faces); s = self.source; n = s.network; sums = amounts.sum(axis=0)
        number = sums[self.number]; mass = s.mass@number/s.neutral_mass; momentum = sums[s.momentum].sum(axis=0)
        particle = sums[s.energy].sum()+n.energy_ev@sums[:n.count]+sums[n.photons]-sums[n.heating]
        particle += (s.momentum_scale/EV)*s.speed_scale*(self.frame@momentum)+.5*mass*(self.frame@self.frame)
        field = self.electric(amounts); scale = EPS0/(6*EV*self.density)
        interior = scale*np.sum(width*(field[:-1]**2+field[:-1]*field[1:]+field[1:]**2))
        exterior = 3*scale*(field[0]**2*(faces[0]-self.bounds[0])+field[-1]**2*(self.bounds[1]-faces[-1]))
        return dict(nuclei=float(number[:-1].sum()), charge=float(self.charge@number), mass_over_neutral=float(mass),
            momentum_scaled=momentum+mass*self.frame, particle_lab_energy_ev=float(particle),
            interior_field_energy_ev=float(interior), exterior_field_energy_ev=float(exterior),
            total_energy_ev=float(particle+interior+exterior))
