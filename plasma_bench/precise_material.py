"""Existing compatible material equations with distinct decimal support geometry."""
from dataclasses import dataclass
from decimal import Decimal
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from scipy.linalg import eigh
from .mass_scaling import solve_mass_banded
from .precise_geometry import coordinates, decimal, widths, fraction, interpolate, geometry_precision
from .species_material import MaterialSpecies, SpeciesMaterialState
from .material_subcells import MaterialSubcells
from .material_species_set import MaterialSpeciesSet


@dataclass
class PreciseMaterialState(SpeciesMaterialState):
    @geometry_precision
    def added(self, rate, scale):
        return PreciseMaterialState(self.faces+decimal(scale)*coordinates(rate.faces),
            self.velocity+scale*rate.velocity, self.thermal+scale*rate.thermal)


class PreciseMaterialSubcells(MaterialSubcells):
    @geometry_precision
    def __init__(self, numbers, mass, charge, bounds, faces, velocity_cells, *, number_scale=None,
                 velocity_faces=None, **controls):
        numbers = np.asarray(numbers, float)
        self.number_scale = float(numbers.max()) if number_scale is None else float(number_scale)
        if not np.isfinite(self.number_scale) or self.number_scale <= 0:
            raise ValueError('Positive finite population coordinate scale required')
        self.relative_numbers = numbers/self.number_scale if number_scale is None else numbers.copy()
        MaterialSpecies.__init__(self, self.relative_numbers[None, :], [mass], [charge],
                                 np.asarray(bounds, float), **controls)
        self.bounds = coordinates(bounds)
        self.numbers = (self.relative_numbers*self.number_scale)[None, :]
        self.reference_faces = coordinates(faces)
        if (self.reference_faces.shape != (self.cells+1,) or np.any(widths(self.reference_faces) <= 0)
                or not isinstance(velocity_cells, int) or velocity_cells < 1):
            raise ValueError('Ordered precise density faces and positive velocity resolution required')
        if velocity_faces is None:
            target = self.reference_faces[0]+(self.reference_faces[-1]-self.reference_faces[0])*coordinates(
                (1-np.cos(np.linspace(0, np.pi, velocity_cells+1)))/2)
            self.velocity_indices = np.unique(np.argmin(abs(self.reference_faces[:, None]-target[None, :]), axis=0))
        else:
            coarse = coordinates(velocity_faces)
            if (coarse.ndim != 1 or len(coarse) < 2 or np.any(widths(coarse) <= 0) or
                    coarse[0] != self.reference_faces[0] or coarse[-1] != self.reference_faces[-1]):
                raise ValueError('Ordered velocity faces spanning the complete occupied support required')
            self.velocity_indices = np.searchsorted(self.reference_faces, coarse)
            if not np.array_equal(self.reference_faces[self.velocity_indices], coarse):
                raise ValueError('Every velocity face must be an existing precise density face')
        self.velocity_cells = len(self.velocity_indices)-1
        coarse = self.reference_faces[self.velocity_indices]
        self.owner = np.searchsorted(coarse, self.reference_faces[:-1], side='right')-1
        length = np.diff(coarse)[self.owner]
        fine = np.stack((self.reference_faces[:-1], self.reference_faces[1:]), axis=1)
        self.exact_right_shape = (fine-coarse[self.owner, None])/length[:, None]
        self.exact_left_shape = Decimal(1)-self.exact_right_shape
        self.span = np.asarray(np.diff(self.reference_faces)/length, float)
        self.left_shape = np.asarray(self.exact_left_shape, float)
        self.right_shape = np.asarray(self.exact_right_shape, float)
        l, r = self.left_shape, self.right_shape
        cell_mass = self.mass[0]*self.relative_numbers
        ll = (l[:, 0]**2+l[:, 0]*l[:, 1]+l[:, 1]**2)/3
        rr = (r[:, 0]**2+r[:, 0]*r[:, 1]+r[:, 1]**2)/3
        lr = (2*l[:, 0]*r[:, 0]+l[:, 0]*r[:, 1]+l[:, 1]*r[:, 0]+2*l[:, 1]*r[:, 1])/6
        matrix = np.zeros((self.velocity_cells+1, self.velocity_cells+1))
        j = self.owner
        for row, col, value in [(j, j, ll), (j+1, j+1, rr), (j, j+1, lr), (j+1, j, lr)]:
            np.add.at(matrix, (row, col), cell_mass*value)
        self.mass_matrix = matrix
        self.mass_scale = float(np.diag(matrix).max())
        self.matrix = np.zeros((1, 3, self.velocity_cells+1))
        self.matrix[0, 1] = np.diag(matrix)
        self.matrix[0, 0, 1:] = np.diag(matrix, 1)
        self.matrix[0, 2, :-1] = np.diag(matrix, -1)

    @geometry_precision
    def validate(self, state):
        if (state.faces.shape != (1, self.cells+1) or state.velocity.shape != (1, self.velocity_cells+1, 3)
                or state.thermal.shape != (1, self.cells) or not np.isfinite(state.velocity).all()
                or not np.isfinite(state.thermal).all() or not all(x.is_finite() for x in state.faces.ravel())):
            raise ValueError('Complete precise geometry and finite physical material state required')
        h = widths(state.faces)
        if np.any(h <= 0) or state.faces[0, 0] <= self.bounds[0] or state.faces[0, -1] >= self.bounds[1]:
            raise ValueError('Precise material cells crossed or reached remote boundary')
        if np.any(state.thermal < 0):
            raise ValueError('Negative precise-subcell thermal energy; no repair')
        return h

    def state(self, faces, temperature_ev, velocity=None):
        velocity = np.zeros((1, self.velocity_cells+1, 3)) if velocity is None else np.asarray(velocity, float).copy()
        result = PreciseMaterialState(coordinates(faces), velocity,
            self.relative_numbers[None, :]*EV*np.broadcast_to(temperature_ev, self.numbers.shape)/(self.gamma-1))
        self.validate(result)
        return result

    @geometry_precision
    def pressure(self, state):
        h = self.validate(state)
        pressure = (self.gamma-1)*state.thermal/h
        sound = np.sqrt(self.gamma*pressure*h/(self.mass[0]*self.relative_numbers[None, :]))
        coarse_h = widths(state.faces[0, self.velocity_indices])
        jump = np.diff(state.velocity[0, :, 2])
        compression = np.maximum(-jump, 0.)
        if self.viscosity_limiter == 'minmod':
            gradient = jump/coarse_h
            left, right = np.r_[gradient[0], gradient[:-1]], np.r_[gradient[1:], gradient[-1]]
            resolved = np.minimum(np.minimum(np.maximum(-left, 0.), np.maximum(-gradient, 0.)), np.maximum(-right, 0.))
            compression = np.maximum(compression-coarse_h*resolved, 0.)
        c = compression[self.owner][None, :]
        q = self.mass[0]*self.relative_numbers[None, :]/h*(self.viscosity[0]*sound*c+self.viscosity[1]*c*c)
        return pressure, q, sound

    @geometry_precision
    def acoustic_timestep(self, state):
        p, q, _ = self.pressure(state)
        h = widths(state.faces[0, self.velocity_indices])
        coefficient = self.gamma*(p[0]+q[0])*self.validate(state)[0]/h[self.owner]**2
        stiffness = np.zeros_like(self.mass_matrix)
        j = self.owner
        for row, col, sign in [(j, j, 1), (j+1, j+1, 1), (j, j+1, -1), (j+1, j, -1)]:
            np.add.at(stiffness, (row, col), sign*coefficient)
        frequency = 0.
        if np.any(coefficient):
            eigen = eigh(stiffness/self.mass_scale, self.mass_matrix/self.mass_scale,
                subset_by_index=[self.velocity_cells, self.velocity_cells], eigvals_only=True, check_finite=False)[0]
            if eigen < 0:
                raise ValueError('Negative maximum acoustic eigenvalue')
            frequency = np.sqrt(eigen)
        deformation = float(np.max(abs(np.diff(state.velocity[0, :, 2]))/h, initial=0.))
        return self.cfl/max(frequency, deformation) if max(frequency, deformation) > 0 else np.inf

    @geometry_precision
    def material_rate(self, state, geometry, frame_z):
        p, q, _ = self.pressure(state)
        stress, j = p+q, self.owner
        force = np.zeros_like(state.velocity)
        weighted = stress[0]*self.span
        np.add.at(force[0, :, 2], j, -weighted)
        np.add.at(force[0, :, 2], j+1, weighted)
        edges, electric, active, _ = geometry
        cell, inside = active[0]
        k = cell[inside]
        owner = self.owner[k]
        coarse = state.faces[0, self.velocity_indices]
        h = widths(edges)[inside]
        nodes = np.array([.5-.5/np.sqrt(3), .5+.5/np.sqrt(3)])
        points = edges[:-1][inside, None]+np.diff(edges)[inside, None]*coordinates(nodes)
        right = fraction(points, coarse[owner, None], coarse[owner+1, None])
        field = electric[:-1][inside, None]+np.diff(electric)[inside, None]*nodes
        weight = EV*self.charge[0]*self.relative_numbers[k]/self.validate(state)[0, k]*h/2
        np.add.at(force[0, :, 2], owner, weight*np.sum(field*(1-right), axis=1))
        np.add.at(force[0, :, 2], owner+1, weight*np.sum(field*right, axis=1))
        acceleration = solve_mass_banded(self.matrix[0], force[0])[None, ...]
        # Keep the small differential motion between fine faces before rounding.
        vz = coordinates(state.velocity[0, :, 2])
        endpoints = self.exact_left_shape*vz[j, None]+self.exact_right_shape*vz[j+1, None]
        face_velocity = np.r_[endpoints[:, 0], endpoints[-1, 1]][None, :]+decimal(frame_z)
        jump = self.span*np.diff(state.velocity[0, :, 2])[j]
        return SpeciesMaterialState(face_velocity, acceleration, -stress*jump[None, :])

    @geometry_precision
    def sample(self, state, points):
        x = coordinates(points)
        h = self.validate(state)[0]
        cell = np.searchsorted(state.faces[0], x, side='right')-1
        inside = (cell >= 0) & (cell < self.cells)
        j = cell[inside]
        owner = self.owner[j]
        coarse = state.faces[0, self.velocity_indices]
        a = fraction(x[inside], coarse[owner], coarse[owner+1])
        number, temperature = np.zeros((1, len(x))), np.zeros((1, len(x)))
        velocity = np.zeros((1, len(x), 3))
        number[0, inside] = self.relative_numbers[j]/h[j]*self.number_scale
        velocity[0, inside] = (1-a[:, None])*state.velocity[0, owner]+a[:, None]*state.velocity[0, owner+1]
        temperature[0, inside] = (self.gamma-1)*state.thermal[0, j]/(self.relative_numbers[j]*EV)
        return dict(number=number, velocity=velocity, temperature_ev=temperature)


class PreciseMaterialSet(MaterialSpeciesSet):
    @geometry_precision
    def field(self, state):
        self.validate(state)
        edges = np.unique(np.r_[self.bounds, state.faces])
        h = widths(edges)
        x = (edges[:-1]+edges[1:])/Decimal(2)
        charge = np.zeros(len(x))
        active = []
        for m, s in zip(self.models, state.blocks):
            j = np.searchsorted(s.faces[0], x, side='right')-1
            inside = (j >= 0) & (j < m.cells)
            active.append((j, inside))
            charge[inside] += EV*m.charge[0]*m.number_density(s)[0, j[inside]]
        electric = np.r_[self.left_electric, self.left_electric+np.cumsum(charge*h)/EPS0]
        energy = EPS0/6*np.sum(h*(electric[:-1]**2+electric[:-1]*electric[1:]+electric[1:]**2))
        return edges, electric, active, float(energy)

    @geometry_precision
    def step(self, state, dt):
        return super().step(state, dt)

    @geometry_precision
    def sample(self, state, points):
        self.validate(state)
        values = [m.sample(s, points) for m, s in zip(self.models, state.blocks)]
        result = {key: np.concatenate([v[key] for v in values], axis=0) for key in ('number', 'velocity', 'temperature_ev')}
        result['velocity'] = np.where(result['number'][..., None] > 0, result['velocity']+self.reference_velocity, 0.)
        result['current'] = EV*np.sum(self.charge[:, None, None]*result['number'][..., None]*result['velocity'], axis=0)
        edges, electric, _, _ = self.field(state)
        result['electric'] = interpolate(points, edges, electric)
        return result

    @geometry_precision
    def totals(self, state):
        self.validate(state)
        number, momentum, kinetic, thermal = [], np.zeros(3), 0., 0.
        for m, s in zip(self.models, state.blocks):
            left, right = (v+self.reference_velocity for v in m.fine_velocity(s.velocity))
            scale, amount = m.number_scale, m.relative_numbers
            mass = m.mass[0]*amount
            number.append(amount.sum()*scale)
            momentum += np.sum(mass[:, None]*(left+right)/2, axis=0)*scale
            kinetic += np.sum(mass/6*np.sum(left*left+left*right+right*right, axis=1))*scale
            thermal += s.thermal.sum()*scale
        edges, electric, active, field = self.field(state)
        empty = ~np.any(np.stack([inside for _, inside in active]), axis=0)
        exterior = EPS0/6*np.sum(widths(edges)[empty]*(electric[:-1][empty]**2+
            electric[:-1][empty]*electric[1:][empty]+electric[1:][empty]**2))
        return dict(numbers_m2=np.asarray(number), momentum_kg_per_m_s=momentum,
            kinetic_j_m2=float(kinetic), thermal_j_m2=float(thermal), electric_energy_j_m2=field,
            empty_field_j_m2=float(exterior), energy_j_m2=float(kinetic+thermal+field))
