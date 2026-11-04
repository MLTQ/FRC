"""Conforming tetrahedral gas mesh for electrostatic experiments (optional Gmsh)."""
from pathlib import Path

import numpy as np

from .iec_vacuum import truncated_cage


def build_mesh(path, *, size=.018, wire_size=.001, cage_diameter=.15,
               wire_diameter=.002, wall_radius=.2, geometry='cage'):
    """Round wires and junctions are excluded volumes, never gas-filled cells."""
    import gmsh
    cage = truncated_cage(cage_diameter, wire_diameter, wall_radius)
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.option.setNumber('General.NumThreads', 2)
        gmsh.model.add('electrostatic')
        occ = gmsh.model.occ
        outer = occ.addSphere(0, 0, 0, wall_radius)
        if geometry == 'cage':
            solids = [(3, occ.addSphere(*p, cage.wire_radius)) for p in cage.vertices]
            for a, b in cage.vertices[cage.edges]:
                solids.append((3, occ.addCylinder(*a, *(b-a), cage.wire_radius)))
            metal, _ = occ.fuse(solids[:1], solids[1:])
            gas, _ = occ.cut([(3, outer)], metal)
        elif geometry == 'shell':
            inner = occ.addSphere(0, 0, 0, cage_diameter/2)
            gas, _ = occ.cut([(3, outer)], [(3, inner)])
        elif geometry == 'sphere':
            gas = [(3, outer)]
        else:
            raise ValueError('Unknown mesh geometry')
        occ.synchronize()
        surfaces = gmsh.model.getBoundary(gas, oriented=False)
        cage_surfaces, wall_surfaces = [], []
        for dim, tag in surfaces:
            bounds = np.asarray(gmsh.model.getBoundingBox(dim, tag))
            # The exterior sphere is the only surface reaching the chamber radius.
            (wall_surfaces if np.max(np.abs(bounds)) > .9*wall_radius else cage_surfaces).append(tag)
        if cage_surfaces:
            distance = gmsh.model.mesh.field.add('Distance')
            gmsh.model.mesh.field.setNumbers(distance, 'SurfacesList', cage_surfaces)
            gmsh.model.mesh.field.setNumber(distance, 'Sampling', 60)
            threshold = gmsh.model.mesh.field.add('Threshold')
            for key, value in dict(InField=distance, SizeMin=wire_size, SizeMax=size,
                                   DistMin=wire_size, DistMax=.025).items():
                gmsh.model.mesh.field.setNumber(threshold, key, value)
            gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
        gmsh.option.setNumber('Mesh.MeshSizeMax', size)
        gmsh.option.setNumber('Mesh.MeshSizeMin', wire_size)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 10)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
        gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
        gmsh.model.mesh.generate(3)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        xyz = xyz.reshape(-1, 3)
        remap = np.full(int(tags.max())+1, -1, dtype=np.int64)
        remap[tags] = np.arange(len(tags))
        _, tet_tags = gmsh.model.mesh.getElementsByType(4)
        tets = remap[tet_tags].reshape(-1, 4)
        triangles, kinds = [], []
        for kind, group in ((1, cage_surfaces), (2, wall_surfaces)):
            for tag in group:
                _, nodes = gmsh.model.mesh.getElementsByType(2, tag)
                tri = remap[nodes].reshape(-1, 3)
                triangles.extend(tri)
                kinds.extend([kind]*len(tri))
        result = dict(nodes=xyz, tets=tets, triangles=np.asarray(triangles, np.int64),
                      boundary_kind=np.asarray(kinds, np.int8), vertices=cage.vertices,
                      edges=cage.edges, wall_radius=np.array(wall_radius),
                      wire_radius=np.array(cage.wire_radius), size=np.array(size),
                      wire_size=np.array(wire_size), geometry=np.array(geometry))
        np.savez_compressed(Path(path), **result)
        return result
    finally:
        gmsh.finalize()
