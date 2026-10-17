"""Conservative, locally bounded nodal moment transfer on retained material cells."""
import numpy as np
from .mass_scaling import solve_mass_banded


def bounded_moments(model, left_momentum, right_momentum, left_energy, right_energy):
    """Synchronize antisymmetric momentum/energy corrections to a positive remap.

    All amounts have the model's existing population scale divided out. This
    changes the numerical transfer, not the physical mass or source moments.
    """
    j = model.owner
    l, r = model.left_shape, model.right_shape
    a, b = l.mean(axis=1), r.mean(axis=1)
    mass = model.mass[0]*model.relative_numbers
    momentum = left_momentum+right_momentum
    energy = left_energy+right_energy
    donor = np.column_stack((momentum/mass[:, None], energy/mass))
    if np.any(donor[:, 3] < .5*np.sum(donor[:, :3]**2, axis=1)):
        raise ValueError('Inadmissible source cell mean; no bounded-transfer repair')
    size = model.velocity_cells+1
    weight = np.zeros(size)
    low = np.zeros((size, 4))
    minimum = np.full((size, 4), np.inf)
    maximum = np.full((size, 4), -np.inf)
    # Positive cell-to-node partition: exact cell N/P/E, including trace cells.
    for node, shape in ((j, a), (j+1, b)):
        share = mass*shape
        np.add.at(weight, node, share)
        np.add.at(low, node, share[:, None]*donor)
        active = shape > 0
        np.minimum.at(minimum, node[active], donor[active])
        np.maximum.at(maximum, node[active], donor[active])
    if np.any(weight <= 0):
        raise ValueError('Unrepresentable nodal mass in bounded reaction transfer')
    # The low state is a convex combination of admissible donor moments.
    specific = low/weight[:, None]
    left = np.column_stack((left_momentum, left_energy))
    right = np.column_stack((right_momentum, right_energy))
    nodal_load = np.zeros((size, 4))
    np.add.at(nodal_load, j, l[:, 0, None]*left+l[:, 1, None]*right)
    np.add.at(nodal_load, j+1, r[:, 0, None]*left+r[:, 1, None]*right)
    high = solve_mass_banded(model.matrix[0], nodal_load)
    flux = model.matrix[0, 0, 1:, None]*(high[:-1]-high[1:])
    # Difference between quadrature and cell-mean loads is itself antisymmetric.
    np.add.at(flux, j, .5*(l[:, 0]-l[:, 1])[:, None]*(left-right))

    def admissible(value):
        return value[3] >= .5*float(value[:3]@value[:3])

    limited = 0
    smallest = 1.
    # Ordered local flux correction. The ordering is declared and measured;
    # each accepted pair conserves all four extensive coordinates.
    for edge, transfer in enumerate(flux):
        changes = (transfer/weight[edge], -transfer/weight[edge+1])
        alpha = 1.
        for node, change in zip((edge, edge+1), changes):
            positive, negative = change > 0, change < 0
            if np.any(positive):
                alpha = min(alpha, float(np.min((maximum[node, positive]-specific[node, positive])/change[positive])))
            if np.any(negative):
                alpha = min(alpha, float(np.min((minimum[node, negative]-specific[node, negative])/change[negative])))
        # A roundoff departure of the low average from its extrema permits no
        # antidiffusion. It does not clip or overwrite a physical state value.
        alpha = max(0., min(1., alpha))
        for node, change in zip((edge, edge+1), changes):
            if not admissible(specific[node]+alpha*change):
                lo, hi = 0., alpha
                for _ in range(54):
                    mid = .5*(lo+hi)
                    if admissible(specific[node]+mid*change):
                        lo = mid
                    else:
                        hi = mid
                alpha = lo
        specific[edge] += alpha*changes[0]
        specific[edge+1] += alpha*changes[1]
        limited += int(alpha < 1.)
        smallest = min(smallest, alpha)
    velocity = specific[:, :3]
    nodal_heat = specific[:, 3]-.5*np.sum(velocity*velocity, axis=1)
    if np.any(nodal_heat < 0):
        raise ValueError('Negative bounded nodal thermal energy; no repair')
    lr = (2*l[:, 0]*r[:, 0]+l[:, 0]*r[:, 1]+l[:, 1]*r[:, 0]+2*l[:, 1]*r[:, 1])/6
    thermal = mass*(a*nodal_heat[j]+b*nodal_heat[j+1]+
        .5*lr*np.sum((velocity[j+1]-velocity[j])**2, axis=1))
    vl, vr = model.fine_velocity(velocity[None, ...])
    kinetic = mass/6*np.sum(vl*vl+vl*vr+vr*vr, axis=1)
    return velocity, thermal, thermal+kinetic-energy, dict(
        limited_edges=limited, total_edges=len(flux), minimum_correction_factor=smallest,
        maximum_nodal_bound_violation=float(max(np.max(minimum-specific), np.max(specific-maximum), 0.)))
