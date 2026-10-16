"""Symmetric reaction/transport coupling of actual variable material moments."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as EV
from .material_reactions import MaterialReactionOverlap
from .material_reaction_projection import project_reactions


@dataclass
class MaterialReservoirs:
    photons_j_m2: float=0.
    input_heat_j_m2: float=0.


class ReactiveMaterial:
    """Source-half / moving-full / source-half, with explicit projection limits."""
    def __init__(self,source,reference_density,*,quadrature=3,**transport_controls):
        self.overlap=MaterialReactionOverlap(source,reference_density,quadrature=quadrature)
        self.controls=transport_controls;self.last_overlap=None

    def source_step(self,model,state,reservoirs,dt):
        gathered=self.overlap.from_material(model,state,model.species_indices)
        self.last_overlap=self.overlap.step(gathered,dt)
        result,advanced,record=project_reactions(self.overlap,self.last_overlap,**self.controls)
        ledger=MaterialReservoirs(reservoirs.photons_j_m2+record['photon_energy_j_m2'],reservoirs.input_heat_j_m2+record['input_heat_j_m2'])
        return result,advanced,ledger

    def step(self,model,state,reservoirs,dt):
        model,state,reservoirs=self.source_step(model,state,reservoirs,dt/2)
        # A reaction can change both pressure and mesh restrictions. Reject a
        # requested interval that no longer resolves the actual transport state.
        if dt>model.timestep(state)*(1+1e-12):raise ValueError('Reaction changed the admissible material timestep; no hidden retry')
        state=model.step(state,dt)
        return self.source_step(model,state,reservoirs,dt/2)

    def totals(self,model,state,reservoirs):
        result=model.totals(state);source=self.overlap.source;number=np.zeros(len(self.overlap.number))
        np.add.at(number,model.species_indices,result['numbers_m2'])
        binding=EV*(source.network.energy_ev@number[:-1])
        result.update(atomic_numbers_m2=number,binding_energy_j_m2=float(binding),photon_energy_j_m2=reservoirs.photons_j_m2,
                      full_energy_j_m2=float(result['energy_j_m2']+binding+reservoirs.photons_j_m2-reservoirs.input_heat_j_m2))
        return result
