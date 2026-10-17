"""Volume-norm velocity projection with physical momentum and cell-energy constraints."""
import numpy as np
from .precise_geometry import fraction
from .mass_scaling import solve_mass_banded


def constrained_velocity(model, points, weights, velocities, occupied, momentum, energy):
    """Minimize volume velocity error without exceeding any cell's total energy.

    A strictly feasible logarithmic-barrier solve enforces the original convex
    constraints. No heat is clipped and no energy reservoir is adjusted.
    """
    coarse = model.reference_faces[model.velocity_indices]
    inside = occupied & (points >= model.reference_faces[0]) & (points <= model.reference_faces[-1])
    x, w, u = points[inside], weights[inside], velocities[inside]
    owner = np.searchsorted(coarse, x, side='right')-1
    owner = np.clip(owner, 0, model.velocity_cells-1)
    right = fraction(x, coarse[owner], coarse[owner+1])
    left = 1-right
    size = model.velocity_cells+1
    matrix = np.zeros((3, size))
    load = np.zeros((size, 3))
    # Normalize only the optimization coordinates, not physical cell states.
    mass = model.mass[0]*model.relative_numbers
    specific_energy = energy/mass
    speed = float(max(np.max(abs(u)), np.sqrt(2*specific_energy.max())))
    if not np.isfinite(speed) or speed <= 0:
        raise ValueError('Positive finite velocity scale required for constrained transfer')
    w = w/w.sum()
    u = u/speed
    np.add.at(matrix[1], owner, w*left*left)
    np.add.at(matrix[1], owner+1, w*right*right)
    np.add.at(matrix[0], owner+1, w*left*right)
    matrix[2, :-1] = matrix[0, 1:]
    np.add.at(load, owner, (w*left)[:, None]*u)
    np.add.at(load, owner+1, (w*right)[:, None]*u)
    physical_weight = model.mass_matrix.sum(axis=1)/mass.sum()
    mean = momentum/mass.sum()/speed
    linear = solve_mass_banded(matrix, np.column_stack((load, physical_weight)))
    v = linear[:, :3]+linear[:, 3, None]*(mean-physical_weight@linear[:, :3])[None, :]/(
        physical_weight@linear[:, 3])
    j = model.owner
    l, r = model.left_shape, model.right_shape
    aa = (l[:, 0]**2+l[:, 0]*l[:, 1]+l[:, 1]**2)/3
    cc = (r[:, 0]**2+r[:, 0]*r[:, 1]+r[:, 1]**2)/3
    bb = (2*l[:, 0]*r[:, 0]+l[:, 0]*r[:, 1]+l[:, 1]*r[:, 0]+2*l[:, 1]*r[:, 1])/6
    ee = specific_energy/speed**2

    def slack(value):
        vl, vr = value[j], value[j+1]
        return ee-.5*(aa*np.sum(vl*vl, axis=1)+2*bb*np.sum(vl*vr, axis=1)+cc*np.sum(vr*vr, axis=1))

    initial_slack = slack(v)
    record = dict(unconstrained_negative_cells=int(np.sum(initial_slack < 0)),
                  barrier_newton_steps=0, barrier_parameter=0.)
    if np.all(initial_slack >= 0):
        return v*speed, record
    # The constant physical mean is a feasible starting point for the present
    # symmetric fixtures. Other states must supply a genuine feasible solution;
    # this routine rejects if that simple starting point is unavailable.
    v = np.broadcast_to(mean, (size, 3)).copy()
    if np.any(slack(v) <= 0):
        raise ValueError('No strict constant-mean feasible start for energy-constrained transfer')
    geom = np.diag(matrix[1])+np.diag(matrix[0, 1:], 1)+np.diag(matrix[2, :-1], -1)
    base_hessian = np.kron(geom, np.eye(3))
    constraint = np.kron(physical_weight[None, :], np.eye(3))
    indices = 3*j[:, None]+np.arange(3)
    next_indices = indices+3
    normalizer = len(j)

    def objective(value, mu):
        s = slack(value)
        if np.any(s <= 0):
            return np.inf
        return .5*np.sum(value*(geom@value))-np.sum(load*value)-mu*np.sum(np.log(s/ee))/normalizer

    for mu in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
        record['barrier_parameter'] = mu
        converged = False
        for _ in range(80):
            s = slack(v)
            dl = aa[:, None]*v[j]+bb[:, None]*v[j+1]
            dr = bb[:, None]*v[j]+cc[:, None]*v[j+1]
            coeff = mu/(normalizer*s)
            gradient = geom@v-load
            np.add.at(gradient, j, coeff[:, None]*dl)
            np.add.at(gradient, j+1, coeff[:, None]*dr)
            hessian = base_hessian.copy()
            for ids1, d1, ids2, d2, q in (
                    (indices, dl, indices, dl, aa), (indices, dl, next_indices, dr, bb),
                    (next_indices, dr, indices, dl, bb), (next_indices, dr, next_indices, dr, cc)):
                local = coeff[:, None, None]*(d1[:, :, None]*d2[:, None, :]/s[:, None, None]+
                    q[:, None, None]*np.eye(3)[None, :, :])
                np.add.at(hessian, (ids1[:, :, None], ids2[:, None, :]), local)
            scale = np.sqrt(np.diag(hessian))
            h = hessian/scale[:, None]/scale[None, :]
            rhs = np.column_stack((gradient.ravel(), constraint.T))/scale[:, None]
            solved = np.linalg.solve(h, rhs)/scale[:, None]
            g, z = solved[:, 0], solved[:, 1:]
            step = -g+z@np.linalg.solve(constraint@z, constraint@g)
            descent = float(gradient.ravel()@step)
            if -.5*descent < max(1e-20, mu*1e-4):
                converged = True
                break
            old = objective(v, mu)
            alpha = 1.
            direction = step.reshape(size, 3)
            for _ in range(60):
                proposal = v+alpha*direction
                if objective(proposal, mu) <= old+1e-4*alpha*descent:
                    v = proposal
                    break
                alpha *= .5
            else:
                raise ValueError('Energy-constrained velocity line search did not converge')
            record['barrier_newton_steps'] += 1
        if not converged:
            raise ValueError('Energy-constrained velocity barrier did not converge')
    record['normalized_momentum_residual'] = float(abs(physical_weight@v-mean).max())
    if np.any(slack(v) <= 0) or record['normalized_momentum_residual'] > 1e-12:
        raise ValueError('Energy-constrained velocity feasibility residual exceeded tolerance')
    return v*speed, record
