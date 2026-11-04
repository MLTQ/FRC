"""Read-only field-line geometry and optional plasma-only magnetic diagnostic."""
import numpy as np
from numba import njit
from plasma_bench.iec_field import walk


def seeds(radius, count):
    k = np.arange(count) + .5
    z = 1-2*k/count
    angle = k*np.pi*(3-np.sqrt(5))
    return radius*np.column_stack((np.sqrt(1-z*z)*np.cos(angle),
                                  np.sqrt(1-z*z)*np.sin(angle), z))


def electric_lines(mesh, electric, spacing=.002, steps=320):
    """Integrate dx/ds = ±E/|E| in native tetrahedra; stop at actual steel."""
    starts = np.concatenate((seeds(.035, 32), seeds(.09, 24), seeds(.155, 16)))
    cells = mesh.locate(starts)
    keep = cells >= 0
    starts, cells = starts[keep], cells[keep]
    result = []
    peak = max(float(np.linalg.norm(electric, axis=1).max()), 1e-30)
    for sign in (-1, 1):
        x, active, cell = starts.copy(), np.ones(len(starts), bool), cells.copy()
        paths = [[p.tolist()] for p in x]
        values = [[float(np.linalg.norm(electric[c]))] for c in cell]
        for step in range(steps):
            ids = np.flatnonzero(active)
            if not len(ids):
                break
            v = electric[cell[ids]]
            speed = np.linalg.norm(v, axis=1)
            moving = speed > peak*1e-10
            active[ids[~moving]] = False
            ids, v, speed = ids[moving], v[moving], speed[moving]
            if not len(ids):
                continue
            middle = x[ids]+sign*.5*spacing*v/speed[:, None]
            mid, midcell, hit = walk(x[ids], middle, cell[ids], mesh.origins,
                                    mesh.inverse, mesh.neighbours, mesh.boundary)
            direction = electric[midcell]
            norm = np.linalg.norm(direction, axis=1)
            target = x[ids]+sign*spacing*direction/np.maximum(norm[:, None], 1e-300)
            target[hit != 0] = mid[hit != 0]
            end, nextcell, endhit = walk(x[ids], target, cell[ids], mesh.origins,
                                       mesh.inverse, mesh.neighbours, mesh.boundary)
            stopped = (hit != 0) | (endhit != 0) | (norm < peak*1e-10)
            # An unresolved traversal is a missing trace, never an invented path.
            bad = (hit == 3) | (endhit == 3)
            for local, identifier in enumerate(ids):
                if not bad[local] and (step % 2 == 0 or stopped[local]):
                    paths[identifier].append(end[local].tolist())
                    values[identifier].append(float(norm[local]))
            x[ids], cell[ids] = end, nextcell
            active[ids[stopped]] = False
        for p, v in zip(paths, values):
            if len(p) > 1:
                if sign < 0:
                    p, v = p[::-1], v[::-1]
                result.append(dict(points=p, values=v))
    return result


@njit(cache=True)
def biot_savart(points, positions, moments, softening):
    """μ0/(4π) Σ (J dV)×R/(R²+ε²)^(3/2), SI; no circuit/metal closure."""
    out = np.zeros((len(points), 3))
    for k in range(len(points)):
        for j in range(len(positions)):
            x, y, z = points[k]-positions[j]
            factor = 1e-7/(x*x+y*y+z*z+softening*softening)**1.5
            a, b, c = moments[j]
            out[k, 0] += factor*(b*z-c*y)
            out[k, 1] += factor*(c*x-a*z)
            out[k, 2] += factor*(a*y-b*x)
    return out


@njit(cache=True)
def grid_value(vector, x, extent):
    n = vector.shape[0]
    q = (x+extent)*(n-1)/(2*extent)
    if np.any(q < 0) or np.any(q >= n-1):
        return np.zeros(3)
    ijk = np.floor(q).astype(np.int64)
    f = q-ijk
    value = np.zeros(3)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                weight = (f[0] if i else 1-f[0])*(f[1] if j else 1-f[1])*(f[2] if k else 1-f[2])
                value += weight*vector[ijk[0]+i, ijk[1]+j, ijk[2]+k]
    return value


def magnetic_lines(nodes, current_moments, softening):
    extent, n = .205, 17
    axis = np.linspace(-extent, extent, n)
    points = np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), axis=-1).reshape(-1, 3)
    nonzero = np.linalg.norm(current_moments, axis=1) > 0
    b = biot_savart(points, nodes[nonzero], current_moments[nonzero], softening).reshape(n, n, n, 3)
    peak = float(np.linalg.norm(b, axis=-1).max())
    paths = []
    for start in np.concatenate((seeds(.04, 16), seeds(.10, 16))):
        for sign in (-1, 1):
            x, path, values = start.copy(), [], []
            for step in range(230):
                v = grid_value(b, x, extent)
                speed = np.linalg.norm(v)
                if speed <= max(peak*1e-8, 1e-30) or np.linalg.norm(x) > .198:
                    break
                if step % 2 == 0:
                    path.append(x.tolist()); values.append(float(speed))
                mid = x+sign*.0015*v/speed
                vmid = grid_value(b, mid, extent)
                x += sign*.003*vmid/max(np.linalg.norm(vmid), 1e-300)
            if len(path) > 1:
                if sign < 0:
                    path, values = path[::-1], values[::-1]
                paths.append(dict(points=path, values=values))
    return dict(lines=paths, max_t=peak, softening_m=softening, grid_points=n,
                scope='Plasma-current contribution only; instantaneous softened Biot–Savart diagnostic. '
                      'No cage/feed/return/displacement currents, steel magnetic response or magnetic force. '
                      'This is not the total magnetic field or an electromagnetic simulation.')
