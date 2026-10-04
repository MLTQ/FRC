"""Compatible skew current/entropy transport on the actual regular point state."""
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.sparse import csr_matrix
from .physical_derivatives import fluid_matrices
from .regular_field_map import periodic
from .variable_compatible_current import weak_vector
from .species_projection import SpeciesProjection


class SpeciesPotential:
    """Density-weighted common scalar trace and its exact Euclidean transpose."""
    def __init__(self, model):
        self.model = model
        self.fluid = model.mapping.fluid
        self.n = model.density_map.point_values
        mesh = model.mesh
        self.radial = .5*(fluid_matrices(self.fluid, mesh.rf, 'left')
                            + fluid_matrices(self.fluid, mesh.rf, 'right'))
        self.axial = csr_matrix(periodic(self.fluid.z.ravel(), mesh.zf[:-1],
                                         self.fluid.length))
        self.denominator = self.sample(self.n)
        if np.any(self.denominator <= 0) or not np.isfinite(self.denominator).all():
            raise ValueError('Positive finite common-trace density required')

    def sample(self, value):
        value = self.model.mapping.validate(value)
        return np.asarray((self.radial @ value.reshape(-1, self.fluid.nz)) @ self.axial.T)

    def apply(self, value):
        return self.sample(self.n*value)/self.denominator

    def adjoint(self, value):
        if np.shape(value) != self.denominator.shape or not np.isfinite(value).all():
            raise ValueError('Finite common-potential covector required')
        return self.n*np.asarray((self.radial.T @ (value/self.denominator))
                                 @ self.axial).reshape(self.n.shape)


class RegularSpecies:
    """Instantaneous bounded pair; frozen density/current transport is skew."""
    def __init__(self, model, *, transfer='projection'):
        self.model = model
        self.mapping = model.mapping
        self.geometry = model.geometry
        self.fluid = self.mapping.fluid
        if transfer not in ('projection', 'point'):
            raise ValueError('Species transfer must be projection or point control')
        self.potential = SpeciesProjection(model) if transfer == 'projection' else SpeciesPotential(model)

    def gradient(self, scalar):
        return self.geometry.sample_electric(
            self.geometry.gradient(self.potential.apply(scalar)))

    def gradient_adjoint(self, vector):
        geo = self.geometry
        er, _, ez = geo.split_electric(weak_vector(geo, vector))
        value = geo.Dx.T @ er.reshape(self.model.mesh.nr, self.model.mesh.nz)
        value += ez.reshape(self.model.mesh.nr+1, self.model.mesh.nz) @ geo.Dz
        return self.potential.adjoint(np.asarray(value))

    def transport(self, current, entropy, *, skew=True):
        """Return K*s (or the preserved unskew candidate), without -a/(e V)."""
        geo = self.geometry
        current = np.asarray(current, float)
        if current.shape != (self.model.mass.shape[0],) or not np.isfinite(current).all():
            raise ValueError('Finite compatible current coefficients required')
        entropy = self.mapping.validate(entropy)
        j = geo.sample_electric(current)
        forward = self.mapping.adjoint(np.sum(j*self.gradient(entropy), axis=-1))
        if not skew:
            return forward
        reverse = self.gradient_adjoint(j*self.mapping.scalar(entropy)[..., None])
        return .5*(forward-reverse)

    def pair(self, state, current):
        state = np.asarray(state, float)
        fields = self.fluid.fields(state)
        if not np.array_equal(state[..., 0], self.model.density_map.point_values):
            raise ValueError('Species and current must share the actual point density')
        current = np.asarray(current, float)
        if current.shape != (self.model.mass.shape[0],) or not np.isfinite(current).all():
            raise ValueError('Finite compatible current coefficients required')
        geo = self.geometry
        rate = np.zeros_like(state)
        forces = []
        potentials = []
        reduced = []
        electric_potentials = []
        for k, (gamma, a) in enumerate(zip(self.fluid.thermal.gammas,
                                         (-self.model.beta, self.model.alpha))):
            s = fields['entropy'][..., k]
            h = fields['specific_internal'][..., k]
            rate[..., k+4] = -a/E*self.transport(current, s)/self.fluid.volume
            # Exact affine separation: constant h or s belongs to the curl kernel.
            ha, sa = float(h.flat[0]), float(s.flat[0])
            dh, ds = h-ha, s-sa
            phi = self.potential.apply(gamma*h-.5*s*h-.5*ha*s+.5*sa*h)
            remainder = -.5*(self.mapping.scalar(dh)[..., None]*self.gradient(ds)
                               - self.mapping.scalar(ds)[..., None]*self.gradient(dh))
            reduced.append(a/E*weak_vector(geo, remainder))
            electric_potentials.append(a/E*geo.gradient(phi))
            forces.append(self.model.mass @ electric_potentials[-1]+reduced[-1])
            potentials.append(phi)
        forcing = sum(forces)
        if not np.isfinite(rate).all() or not np.isfinite(forcing).all():
            raise ValueError('Species rate exceeds finite numerical range')
        return dict(fluid_rate=rate, forcing=forcing,
                    reduced_forcing=sum(reduced), potential_electric=sum(electric_potentials),
                    species_forcing=np.stack(forces), potentials=np.stack(potentials))
