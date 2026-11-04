"""Linear tetrahedral Poisson solve and charge/field interpolation."""
import numpy as np
from scipy.constants import epsilon_0
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree
from numba import njit


@njit(cache=True)
def barycentric(point, origin, inverse):
    tail = inverse @ (point-origin)
    return np.array((1-tail.sum(), tail[0], tail[1], tail[2]))


@njit(cache=True)
def walk(positions, targets, cells, origins, inverse, neighbours, boundary):
    """Trace straight drifts through adjacent tetrahedra to first metal hit."""
    result = cells.copy()
    hit = np.zeros(len(cells), np.int8)
    out = targets.copy()
    for p in range(len(cells)):
        cell = cells[p]
        start = positions[p].copy()
        for _ in range(1000):
            end_b = barycentric(targets[p], origins[cell], inverse[cell])
            if end_b.min() >= -1e-10:
                result[p] = cell
                break
            begin_b = barycentric(start, origins[cell], inverse[cell])
            fraction = 2.
            face = -1
            for j in range(4):
                if end_b[j] < -1e-10:
                    f = max(0., begin_b[j])/(max(0., begin_b[j])-end_b[j])
                    if f < fraction:
                        fraction, face = f, j
            if face < 0:
                hit[p] = 3
                break
            start += fraction*(targets[p]-start)
            nxt = neighbours[cell, face]
            if nxt < 0:
                hit[p] = boundary[cell, face]
                out[p] = start
                result[p] = cell
                break
            cell = nxt
        else:
            hit[p] = 3
    return out, result, hit


class Field:
    def __init__(self, mesh):
        self.nodes, self.tets = mesh['nodes'], mesh['tets']
        x = self.nodes[self.tets]
        self.origins = x[:, 0].copy()
        jacobian = np.transpose(x[:, 1:]-x[:, :1], (0, 2, 1))
        self.inverse = np.linalg.inv(jacobian)
        self.volume = np.abs(np.linalg.det(jacobian))/6
        self.gradient = np.concatenate((-self.inverse.sum(axis=1)[:, None], self.inverse), axis=1)
        local = epsilon_0*self.volume[:, None, None]*np.einsum('tix,tjx->tij', self.gradient, self.gradient)
        n = len(self.nodes)
        self.K = coo_matrix((local.ravel(), (np.repeat(self.tets, 4, axis=1).ravel(),
                            np.tile(self.tets, (1, 4)).ravel())), shape=(n, n)).tocsr()
        self.kind = np.zeros(n, np.int8)
        for tri, kind in zip(mesh['triangles'], mesh['boundary_kind']):
            self.kind[tri] = kind
        self.free = np.flatnonzero(self.kind == 0)
        self.fixed = np.flatnonzero(self.kind != 0)
        self.factor = splu(self.K[self.free][:, self.free].tocsc())
        self.unit = np.zeros(n)
        self.unit[self.kind == 1] = 1.
        self.unit[self.free] = self.factor.solve(-(self.K@self.unit)[self.free])
        self.capacitance = float(self.unit@(self.K@self.unit))
        self.tree = cKDTree(x.mean(axis=1))
        self.neighbours = np.full((len(self.tets), 4), -1, np.int64)
        self.boundary = np.zeros_like(self.neighbours, dtype=np.int8)
        faces = np.array([np.delete(np.arange(4), j) for j in range(4)])
        keys = np.sort(self.tets[:, faces].reshape(-1, 3), axis=1)
        order = np.lexsort(keys.T[::-1])
        sorted_keys = keys[order]
        pairs = np.flatnonzero(np.all(sorted_keys[:-1] == sorted_keys[1:], axis=1))
        left, right = order[pairs], order[pairs+1]
        self.neighbours.ravel()[left] = right//4
        self.neighbours.ravel()[right] = left//4
        boundary_map = {tuple(sorted(tri)): int(kind) for tri, kind in zip(mesh['triangles'], mesh['boundary_kind'])}
        for flat in np.flatnonzero(self.neighbours.ravel() < 0):
            self.boundary.ravel()[flat] = boundary_map[tuple(keys[flat])]
        self.mesh = mesh

    def locate(self, points):
        # Bounded exhaustive fallback is only for initial loading/slice sampling.
        _, near = self.tree.query(points, k=min(24, len(self.tets)))
        near = np.asarray(near).reshape(len(points), -1)
        tail = np.einsum('pkij,pkj->pki', self.inverse[near], points[:, None]-self.origins[near])
        weights = np.concatenate((1-tail.sum(axis=2, keepdims=True), tail), axis=2)
        inside = weights.min(axis=2) >= -1e-9
        cells = np.where(inside.any(axis=1), near[np.arange(len(points)), inside.argmax(axis=1)], -1)
        for p in np.flatnonzero(cells < 0):
            tail = np.einsum('tij,tj->ti', self.inverse, points[p]-self.origins)
            ok = (tail.min(axis=1) >= -1e-9) & (tail.sum(axis=1) <= 1+1e-9)
            if ok.any():
                cells[p] = np.flatnonzero(ok)[0]
        return cells

    def weights(self, positions, cells):
        tail = np.einsum('pij,pj->pi', self.inverse[cells], positions-self.origins[cells])
        return np.concatenate((1-tail.sum(axis=1, keepdims=True), tail), axis=1)

    def deposit(self, positions, cells, charge):
        return np.bincount(self.tets[cells].ravel(),
                           (charge[:, None]*self.weights(positions, cells)).ravel(), minlength=len(self.nodes))

    def solve(self, charge, voltage):
        phi_sc = np.zeros(len(self.nodes))
        phi_sc[self.free] = self.factor.solve(charge[self.free])
        phi = phi_sc+voltage*self.unit
        electric = -np.einsum('ti,tix->tx', phi[self.tets], self.gradient)
        metal_charge = self.K@phi-charge
        electrodes = np.array([metal_charge[self.kind == k].sum() for k in (1, 2)])
        return phi, electric, electrodes, .5*float(phi@(self.K@phi)), .5*float(charge@phi_sc)

    def drift(self, positions, targets, cells):
        return walk(positions, targets, cells, self.origins, self.inverse, self.neighbours, self.boundary)
