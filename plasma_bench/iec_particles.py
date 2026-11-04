"""Kinetic electron/ion population and physical particle-boundary sources."""
from dataclasses import dataclass

import numpy as np
from scipy.constants import elementary_charge as qe, electron_mass, proton_mass


@dataclass
class Particles:
    x: np.ndarray
    v: np.ndarray
    cells: np.ndarray
    species: np.ndarray
    weight: float
    ion_mass: float

    @property
    def charge(self):
        return np.where(self.species == 0, -qe, qe)*self.weight

    @property
    def mass(self):
        return np.where(self.species == 0, electron_mass, self.ion_mass)*self.weight

    def kinetic(self):
        return .5*self.mass*np.sum(self.v*self.v, axis=1)

    def retain(self, keep):
        for name in ('x', 'v', 'cells', 'species'):
            setattr(self, name, getattr(self, name)[keep])


def load_plasma(field, rng, count, radius, density, te, ti, ion_amu=2.):
    """Specified neutral initial plasma, independent random 3D species positions."""
    positions = []
    for _ in range(2):
        accepted = []
        total = 0
        while total < count:
            x = rng.uniform(-radius, radius, (max(count, 64), 3))
            x = x[np.linalg.norm(x, axis=1) < radius]
            cells = field.locate(x)
            x = x[cells >= 0]
            accepted.append(x)
            total += len(x)
        positions.append(np.concatenate(accepted)[:count])
    x = np.concatenate(positions)
    species = np.repeat(np.arange(2, dtype=np.int8), count)
    ion_mass = ion_amu*proton_mass
    thermal = np.sqrt(qe*np.where(species == 0, te/electron_mass, ti/ion_mass))
    v = rng.normal(size=x.shape)*thermal[:, None]
    return Particles(x, v, field.locate(x), species, density*4*np.pi*radius**3/(3*count), ion_mass)


class WallEmitter:
    """Idealized uniform grounded-wall electron emitter, not gas ionization."""
    def __init__(self, field):
        wall = []
        cells = []
        for cell, face in zip(*np.where(field.boundary == 2)):
            wall.append(np.delete(field.tets[cell], face))
            cells.append(cell)
        self.tri = field.nodes[np.asarray(wall)]
        self.cells = np.asarray(cells)
        cross = np.cross(self.tri[:, 1]-self.tri[:, 0], self.tri[:, 2]-self.tri[:, 0])
        area = np.linalg.norm(cross, axis=1)/2
        self.probability = area/area.sum()
        self.normals = cross/(2*area[:, None])
        centres = field.nodes[field.tets[self.cells]].mean(axis=1)
        inward = centres-self.tri.mean(axis=1)
        self.normals *= np.where(np.sum(self.normals*inward, axis=1) > 0, 1., -1.)[:, None]
        self.centres = centres

    def inject(self, particles, rng, current, dt, energy_ev):
        count = rng.poisson(current*dt/(qe*particles.weight))
        if count == 0:
            return 0., 0.
        select = rng.choice(len(self.tri), count, p=self.probability)
        weights = rng.exponential(size=(count, 3))
        weights /= weights.sum(axis=1)[:, None]
        x = np.einsum('pi,pij->pj', weights, self.tri[select])
        # Roundoff-only displacement into the face's gas tetra (not a sheath).
        x += 1e-9*(self.centres[select]-x)
        v = np.sqrt(2*qe*energy_ev/electron_mass)*self.normals[select]
        for name, values in (('x', x), ('v', v), ('cells', self.cells[select]),
                             ('species', np.zeros(count, np.int8))):
            setattr(particles, name, np.concatenate((getattr(particles, name), values)))
        return count*particles.weight*qe*energy_ev, count*particles.weight*qe
