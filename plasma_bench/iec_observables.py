"""Saved-state full-population observables for the local IEC/HOPS viewer."""
from collections import OrderedDict
from pathlib import Path
import json
import re
import threading

import numpy as np
from scipy.constants import elementary_charge as qe
from scipy.spatial import cKDTree

from plasma_bench.iec_view_fields import electric_lines, magnetic_lines


class ViewMesh:
    """Geometry/interpolation only: no Poisson factorization or field evolution."""
    def __init__(self, path):
        with np.load(path) as a:
            self.nodes, self.tets = a['nodes'].copy(), a['tets'].copy()
            self.size = float(a['size'])
            triangles, kinds = a['triangles'].copy(), a['boundary_kind'].copy()
        x = self.nodes[self.tets]
        self.origins = x[:, 0].copy()
        jacobian = np.transpose(x[:, 1:]-x[:, :1], (0, 2, 1))
        self.inverse = np.linalg.inv(jacobian)
        volume = np.abs(np.linalg.det(jacobian))/6
        self.gradient = np.concatenate((-self.inverse.sum(axis=1)[:, None], self.inverse), axis=1)
        self.dual_volume = np.bincount(self.tets.ravel(), np.repeat(volume/4, 4), minlength=len(self.nodes))
        self.tree = cKDTree(x.mean(axis=1))
        self.neighbours = np.full(self.tets.shape, -1, np.int64)
        self.boundary = np.zeros(self.tets.shape, np.int64)
        faces = np.array([np.delete(np.arange(4), j) for j in range(4)])
        keys = np.sort(self.tets[:, faces].reshape(-1, 3), axis=1)
        order = np.lexsort(keys.T[::-1]); sorted_keys = keys[order]
        pairs = np.flatnonzero(np.all(sorted_keys[:-1] == sorted_keys[1:], axis=1))
        left, right = order[pairs], order[pairs+1]
        self.neighbours.ravel()[left], self.neighbours.ravel()[right] = right//4, left//4
        boundary_map = {tuple(sorted(t)): int(k) for t, k in zip(triangles, kinds)}
        for flat in np.flatnonzero(self.neighbours.ravel() < 0):
            self.boundary.ravel()[flat] = boundary_map[tuple(keys[flat])]
        self.axis = np.linspace(-.18, .18, 41)
        gx, gz = np.meshgrid(self.axis, self.axis)
        self.slice_points = np.column_stack((gx.ravel(), np.zeros(gx.size), gz.ravel()))
        self.slice_cells = self.locate(self.slice_points)
        self.slice_valid = self.slice_cells >= 0
        self.slice_weights = self.weights(self.slice_points[self.slice_valid], self.slice_cells[self.slice_valid])
        g = np.linspace(-.17, .17, 23)
        grid = np.stack(np.meshgrid(g, g, g, indexing='ij'), axis=-1).reshape(-1, 3)
        grid = grid[np.linalg.norm(grid, axis=1) < .18]
        cells = self.locate(grid); keep = cells >= 0
        self.grid, self.grid_cells = grid[keep], cells[keep]
        self.grid_weights = self.weights(self.grid, self.grid_cells)

    def weights(self, points, cells):
        tail = np.einsum('pij,pj->pi', self.inverse[cells], points-self.origins[cells])
        return np.column_stack((1-tail.sum(axis=1), tail))

    def locate(self, points):
        _, near = self.tree.query(points, k=24)
        tail = np.einsum('pkij,pkj->pki', self.inverse[near], points[:, None]-self.origins[near])
        inside = (tail.min(axis=2) >= -1e-9) & (tail.sum(axis=2) <= 1+1e-9)
        cell = np.where(inside.any(axis=1), near[np.arange(len(points)), inside.argmax(axis=1)], -1)
        for k in np.flatnonzero(cell < 0):
            if np.linalg.norm(points[k]) >= .2:
                continue
            t = np.einsum('tij,tj->ti', self.inverse, points[k]-self.origins)
            ok = (t.min(axis=1) >= -1e-9) & (t.sum(axis=1) <= 1+1e-9)
            if ok.any(): cell[k] = np.flatnonzero(ok)[0]
        return cell

    def moment(self, x, cells, values):
        values = np.asarray(values)
        return np.bincount(self.tets[cells].ravel(),
                           (values[:, None]*self.weights(x, cells)).ravel(), minlength=len(self.nodes))

    def interpolate(self, nodal, cells, weights):
        return np.einsum('pi,pi...->p...', weights, nodal[self.tets[cells]])


def moments(mesh, x, v, species, cells, weight):
    density, current = [], []
    for s in (0, 1):
        select = species == s
        number = mesh.moment(x[select], cells[select], np.full(np.sum(select), weight))
        flux = np.column_stack([mesh.moment(x[select], cells[select], weight*v[select, k]) for k in range(3)])
        density.append(number/mesh.dual_volume)
        current.append((2*s-1)*qe*flux/mesh.dual_volume[:, None])
    return density, current


def snapshot(mesh, path, magnetic=False):
    with np.load(path) as a:
        phi, x, v, species, cells = (a[k].copy() for k in ('phi','x','v','species','cells'))
        weight, time_us = float(a['macro_weight']), float(a['time_s'])*1e6
    electric = -np.einsum('ti,tix->tx', phi[mesh.tets], mesh.gradient)
    density, current = moments(mesh, x, v, species, cells, weight)
    je, ji = current; ne, ni = density
    nodal = dict(ne=ne, ni=ni, rho=qe*(ni-ne), je=np.linalg.norm(je, axis=1),
                 ji=np.linalg.norm(ji, axis=1), j=np.linalg.norm(je+ji, axis=1))
    slices = {}
    for key, values in nodal.items():
        a = np.full(len(mesh.slice_points), np.nan)
        a[mesh.slice_valid] = mesh.interpolate(values, mesh.slice_cells[mesh.slice_valid], mesh.slice_weights)
        slices[key] = [None if not np.isfinite(t) else float(t) for t in a]
    e = np.full(len(mesh.slice_points), np.nan)
    e[mesh.slice_valid] = np.linalg.norm(electric[mesh.slice_cells[mesh.slice_valid]], axis=1)
    slices['electric'] = [None if not np.isfinite(t) else float(t) for t in e]
    grid_current = mesh.interpolate(je+ji, mesh.grid_cells, mesh.grid_weights)
    speed = np.linalg.norm(grid_current, axis=1)
    active = speed > max(float(speed.max())*1e-4, 1e-30)
    ion = species == 1
    if ion.any():
        dx, dv = x[ion]-x[ion].mean(axis=0), v[ion]-v[ion].mean(axis=0)
        radius = np.linalg.norm(dx, axis=1)
        vr = np.sum(dx*dv, axis=1)/np.maximum(radius, 1e-30)
        phase = np.column_stack((radius, vr)).tolist()
    else:
        phase = []
    edges = np.linspace(0, .20, 25)
    shell = 4*np.pi/3*np.diff(edges**3)
    profiles = [np.histogram(np.linalg.norm(x[species == s], axis=1), edges)[0]*weight/shell for s in (0,1)]
    integral = [np.sum(j*mesh.dual_volume[:, None], axis=0) for j in current]
    direct = [np.sum(v[species==s],axis=0)*weight*(2*s-1)*qe for s in (0,1)]
    result = dict(time_us=time_us, source='Full native state; display particle count does not enter',
                  electric=dict(lines=electric_lines(mesh,electric), max_vm=float(np.linalg.norm(electric,axis=1).max())),
                  current=dict(points=mesh.grid[active].tolist(), vectors=grid_current[active].tolist(),
                               max_am2=float(speed.max()),scope='Instantaneous particle velocity moment; not a charge-conserving trajectory current.'),
                  slice=dict(axis_m=mesh.axis.tolist(),values=slices), phase=phase,
                  radial=dict(radius_m=((edges[:-1]+edges[1:])/2).tolist(),ne=profiles[0].tolist(),ni=profiles[1].tolist(),
                              scope='Spherical shell diagnostic about chamber origin; no symmetry imposed on particles.'),
                  integrals=dict(number=[float(np.sum(n*mesh.dual_volume)) for n in density],
                                 particle_number=[int(np.sum(species==s))*weight for s in (0,1)],
                                 current_am=[j.tolist() for j in integral],particle_current_am=[j.tolist() for j in direct]))
    if magnetic:
        result['magnetic'] = magnetic_lines(mesh.nodes, (je+ji)*mesh.dual_volume[:,None], mesh.size)
    return result


class ObservableStore:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.meshes, self.frames = OrderedDict(), OrderedDict()
        self.lock = threading.Lock()

    def get(self, identifier, index, magnetic=False):
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,100}',identifier):
            raise ValueError('Invalid run identifier')
        run = self.root/'outputs/lab_runs'/identifier
        if not (run/'run.json').is_file():
            run = self.root/'outputs'/identifier
        run = run.resolve()
        if not run.is_relative_to(self.root/'outputs'):
            raise ValueError('Run must be inside local scientific outputs')
        if not (run/'run.json').is_file() or not (run/'mesh.npz').is_file():
            raise ValueError('No native electrostatic run')
        manifest = json.loads((run/'run.json').read_text())
        if not 0 <= index < manifest.get('frames',0):
            raise ValueError('No saved state at this index')
        state = run/f'state_{index:04d}.npz'
        key = (str(state), state.stat().st_mtime_ns, magnetic)
        with self.lock:
            if key in self.frames:
                self.frames.move_to_end(key); return self.frames[key]
            meshpath = run/'mesh.npz'
            # Mesh objects remain read-only; bounded memory and serialized work.
            meshkey = (str(meshpath),meshpath.stat().st_mtime_ns)
            if meshkey not in self.meshes:
                self.meshes[meshkey] = ViewMesh(meshpath)
                while len(self.meshes)>2: self.meshes.popitem(last=False)
            result = snapshot(self.meshes[meshkey],state,magnetic)
            self.frames[key] = result
            while len(self.frames)>8: self.frames.popitem(last=False)
            return result
