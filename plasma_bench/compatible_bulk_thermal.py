"""Bulk log-entropy advection through the actual conservative continuity flux."""
from dataclasses import dataclass
import numpy as np
from .compatible_particle_flux import CompatibleParticleFlux,ParticleFaces


@dataclass(frozen=True)
class BulkThermalRate:
    entropy: np.ndarray
    density: np.ndarray
    face_work: ParticleFaces
    residuals: dict


class CompatibleBulkThermal:
    """Thermal and face-work component; momentum still belongs to its caller."""
    def __init__(self,thermal,particle_flux=None):
        self.thermal=thermal;self.scalar=thermal.scalar;self.geometry=thermal.geometry
        self.flux=particle_flux if particle_flux is not None else CompatibleParticleFlux(thermal.model.density_map)
        if self.flux.density_map is not thermal.model.density_map:
            raise ValueError('Bulk thermal and continuity must share the actual density map')

    def rate(self,entropy,faces):
        scalar=self.scalar;thermal=self.thermal
        fields=thermal.fields(entropy);gradient=scalar.gradient(entropy)
        particle=self.flux.sample(faces);nt=self.flux.density_rate(faces)
        st=scalar.solve(-scalar.adjoint(np.sum(particle*gradient,axis=-1)))
        h=scalar.project(fields['internal']/thermal.model.density)
        advection=self.flux.vector_adjoint(scalar.sample(h.values)[...,None]*gradient)
        compression=self.flux.density_adjoint(thermal.density_energy_gradient(entropy))
        work=ParticleFaces(compression.radial-advection.radial,compression.axial-advection.axial)
        return BulkThermalRate(st.values,nt,work,dict(entropy_mass=st.residual,
            thermal_projection=h.residual,entropy_iterations=st.iterations,projection_iterations=h.iterations))
