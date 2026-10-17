"""Conservative velocity-shape transfer with explicit nonnegative mixing heat."""
import numpy as np
from .mass_scaling import solve_mass_banded
from .precise_geometry import fraction
from .precise_search import searchsorted


def variance_projection(model, points, weights, source_velocity, relative_number,
                        source_heat_density, groups, momentum, energy):
    """Pair a volume-fitted velocity shape with directly integrated source heat."""
    inside = ((relative_number > 0) & (points >= model.reference_faces[0]) &
              (points <= model.reference_faces[-1]))
    x, weights, u = points[inside], weights[inside], source_velocity[inside]
    number, heat, cell = relative_number[inside], source_heat_density[inside], groups[inside]
    if np.any(heat < 0):
        raise ValueError('Negative source thermal energy in variance projection; no repair')
    coarse = model.reference_faces[model.velocity_indices]
    owner = np.clip(searchsorted(coarse, x, side='right')-1, 0, model.velocity_cells-1)
    right = fraction(x, coarse[owner], coarse[owner+1])
    left = 1-right
    mass, total_number = model.mass[0], model.relative_numbers.sum()
    mean = momentum/(mass*total_number)
    source_mass_weights = weights*number/total_number
    target_mass_weights = model.relative_numbers/total_number
    nodal_mass_weights = model.mass_matrix.sum(axis=1)/(mass*total_number)
    speed = float(np.max(abs(u-mean)))
    size = model.velocity_cells+1
    factor, source_norm, target_norm, candidate_factor = 0., 0., 0., 0.
    if speed == 0:
        velocity = np.broadcast_to(mean, (size, 3)).copy()
    else:
        w = weights/weights.sum()
        matrix, load = np.zeros((3, size)), np.zeros((size, 3))
        centered = (u-mean)/speed
        np.add.at(matrix[1], owner, w*left*left)
        np.add.at(matrix[1], owner+1, w*right*right)
        np.add.at(matrix[0], owner+1, w*left*right)
        matrix[2, :-1] = matrix[0, 1:]
        np.add.at(load, owner, (w*left)[:, None]*centered)
        np.add.at(load, owner+1, (w*right)[:, None]*centered)
        solved = solve_mass_banded(matrix, np.column_stack((load, nodal_mass_weights)))
        shape = solved[:, :3]-solved[:, 3, None]*(nodal_mass_weights@solved[:, :3])[None, :]/(
            nodal_mass_weights@solved[:, 3])
        shape -= nodal_mass_weights@shape
        sampled = left[:, None]*shape[owner]+right[:, None]*shape[owner+1]
        vl, vr = model.fine_velocity(shape[None, ...])
        source_norm = float(np.sum(source_mass_weights*np.sum(sampled*sampled, axis=1)))
        target_norm = float(np.sum(target_mass_weights*np.sum(vl*vl+vl*vr+vr*vr, axis=1))/3)
        denominator = source_norm+target_norm
        if denominator > 0:
            candidate_factor = 2*float(np.sum(source_mass_weights*np.sum(sampled*centered, axis=1)))/denominator
            # Both zero and this second root conserve energy. Choose the root
            # with smaller occupied-volume velocity error, not a distant ray root.
            constant_error = np.sum(w*np.sum(centered**2, axis=1))
            candidate_error = np.sum(w*np.sum((centered-candidate_factor*sampled)**2, axis=1))
            if candidate_error < constant_error:
                factor = candidate_factor
        elif denominator != 0:
            raise ValueError('Invalid velocity-shape norm')
        velocity = mean+speed*factor*shape
    if not np.isfinite(velocity).all() or not np.isfinite(factor):
        raise ValueError('Unrepresentable conservative velocity shape')
    sampled_velocity = left[:, None]*velocity[owner]+right[:, None]*velocity[owner+1]
    mixing_density = .5*mass*number*np.sum((u-sampled_velocity)**2, axis=1)
    thermal = np.zeros(model.cells)
    np.add.at(thermal, cell, weights*(heat+mixing_density))
    vl, vr = model.fine_velocity(velocity[None, ...])
    kinetic = mass*model.relative_numbers/6*np.sum(vl*vl+vl*vr+vr*vr, axis=1)
    redistribution = thermal+kinetic-energy
    energy_scale = float(np.sum(abs(energy)))
    energy_residual = (float(redistribution.sum()/energy_scale) if energy_scale > 0
                       else (0. if np.all(redistribution == 0) else np.inf))
    record = dict(velocity_shape_factor=factor, candidate_shape_factor=candidate_factor, source_shape_norm=source_norm,
        target_shape_norm=target_norm,
        normalized_momentum_residual=float(abs(nodal_mass_weights@velocity-mean).max()/
            max(float(abs(velocity).max()), np.finfo(float).tiny)),
        relative_block_energy_change=energy_residual)
    if (not np.isfinite(thermal).all() or np.any(thermal < 0) or
            record['normalized_momentum_residual'] > 1e-12 or
            abs(record['relative_block_energy_change']) > 1e-12):
        raise ValueError('Conservative variance transfer exceeded physical residual tolerance')
    return velocity, thermal, redistribution, record
