"""Metal transport provider for the shared conservative electron-energy update."""
from .hybrid3d_energy_electrons import EnergyElectrons
from .hybrid3d_metal import MetalKernels


class MetalEnergyElectrons(EnergyElectrons):
    def __init__(self, model, **options):
        super().__init__(model, **options)
        self.metal = MetalKernels(model)
        self.metal.prepare_geometry(self)

    def transport_fluxes(self, u, velocity):
        return self.metal.transport(u, velocity)
