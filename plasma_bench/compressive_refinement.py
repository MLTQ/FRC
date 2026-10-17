"""Bisect full material cells at explicitly unresolved compression."""
import numpy as np
from scipy.constants import elementary_charge as EV
from .species_material import MaterialSpecies, SpeciesMaterialState
from .material_species_set import MaterialSpeciesSet, MaterialSetState


def refine_compressive_cells(model, state, threshold, minimum_child_width, *, maximum_cells=8192, species=None,
                             thermal_reconstruction='constant'):
    if (not np.isfinite([threshold, minimum_child_width]).all() or threshold <= 0 or minimum_child_width <= 0
            or not isinstance(maximum_cells, int) or maximum_cells < 2):
        raise ValueError('Explicit positive refinement threshold, width and integer capacity required')
    if thermal_reconstruction not in ('constant', 'limited_linear'):
        raise ValueError('Unknown thermal reconstruction')
    selected = set(model.species_indices if species is None else species)
    models, blocks, records = [], [], []
    changed = False
    for index, component, block in zip(model.species_indices, model.models, state.blocks):
        if type(component) is not MaterialSpecies or component.species != 1:
            raise ValueError('This refinement requires full-cell single-species material blocks')
        width = component.validate(block)[0]
        jump = np.diff(block.velocity[0, :, 2])
        gradient = jump/width
        compression = np.maximum(-gradient, 0.)
        resolved = np.minimum(np.minimum(np.r_[compression[0], compression[:-1]], compression),
                              np.r_[compression[1:], compression[-1]])
        residual = np.maximum(0., np.maximum(-jump, 0.)-width*resolved)
        sound = component.pressure(block)[2][0]
        requested = (residual > threshold*sound) & (index in selected)
        split = requested & (width/2 >= minimum_child_width)
        added = int(split.sum())
        records.append(dict(species=int(index), initial_cells=component.cells, added_cells=added,
                            width_limited_cells=int(np.count_nonzero(requested & ~split))))
        if not added:
            models.append(component); blocks.append(block); continue
        if component.cells+added > maximum_cells:
            raise ValueError('Requested compression refinement exceeds explicit cell capacity')
        count = 1+split.astype(int)
        offsets = np.r_[0, np.cumsum(count)]
        faces = np.empty(offsets[-1]+1)
        velocity = np.empty((len(faces), 3))
        faces[offsets] = block.faces[0]
        velocity[offsets] = block.velocity[0]
        middle = offsets[:-1][split]+1
        faces[middle] = block.faces[0, :-1][split]+width[split]/2
        velocity[middle] = (block.velocity[0, :-1][split]+block.velocity[0, 1:][split])/2
        numbers = np.repeat(component.numbers[0]/count, count)
        thermal = np.repeat(block.thermal[0]/count, count)
        correction = np.zeros(component.cells)
        if thermal_reconstruction == 'limited_linear' and component.cells > 2:
            pressure = (component.gamma-1)*block.thermal[0]/width
            centre = block.faces[0, :-1]+width/2
            secant = np.diff(pressure)/np.diff(centre)
            slope = np.zeros(component.cells)
            left, right = secant[:-1], secant[1:]
            same_sign = (np.sign(left) == np.sign(right))
            slope[1:-1] = np.where(same_sign, np.sign(left)*np.minimum(abs(left), abs(right)), 0.)
            correction = slope*width**2/(8*(component.gamma-1))
            thermal[offsets[:-1][split]] -= correction[split]
            thermal[offsets[:-1][split]+1] += correction[split]
        records[-1].update(
            thermal_redistribution_l1_j_m2=float(2*abs(correction[split]).sum()),
            maximum_temperature_change_ev=float(np.max(
                (component.gamma-1)*abs(correction[split])/numbers[offsets[:-1][split]]/EV)))
        new_model = MaterialSpecies(numbers[None, :], component.mass, component.charge, component.bounds,
            left_electric=component.left_electric, cfl=component.cfl,
            maximum_plasma_phase=component.maximum_plasma_phase,
            numerical_viscosity=component.viscosity, viscosity_limiter=component.viscosity_limiter)
        new_state = SpeciesMaterialState(faces[None, :], velocity[None, :], thermal[None, :])
        new_model.validate(new_state)
        models.append(new_model); blocks.append(new_state); changed = True
    record = dict(blocks=records, threshold=threshold, minimum_child_width_m=minimum_child_width,
                  maximum_cells=maximum_cells, changed=changed, thermal_reconstruction=thermal_reconstruction)
    if not changed:
        return model, state, record
    result = MaterialSpeciesSet(models, model.species_indices, reference_velocity=model.reference_velocity)
    updated = MaterialSetState(tuple(blocks))
    before, after = model.totals(state), result.totals(updated)
    record.update(relative_energy_change=(after['energy_j_m2']-before['energy_j_m2'])/before['energy_j_m2'],
                  maximum_relative_number_change=float(np.max(abs(after['numbers_m2']-before['numbers_m2'])/before['numbers_m2'])),
                  momentum_change_kg_per_m_s=(after['momentum_kg_per_m_s']-before['momentum_kg_per_m_s']).tolist())
    return result, updated, record
