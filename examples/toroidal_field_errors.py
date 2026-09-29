"""Physical-volume errors of actual reconstruction against continuum point fields."""
import numpy as np
from plasma_bench.toroidal_state import primitive
from plasma_bench.toroidal_reconstruction import ToroidalReconstruction
from .toroidal_references import point_fields


def physical_errors(case, solver, order):
    mesh = solver.mesh
    q = primitive(mesh, solver.state, case.gamma)
    reconstruction = ToroidalReconstruction(mesh, q, order=case.order)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    numerator, denominator = np.zeros(3), np.zeros(3)
    for xr, wr in zip(nodes, weights):
        radius = mesh.r+.5*mesh.dr*xr
        r = radius[:, None]
        for xz, wz in zip(nodes, weights):
            z = mesh.z[None, :]+.5*mesh.dz*xz
            values = reconstruction.raw(radius, .5*xz)
            measured = (values[..., 0], r*values[..., 6], r*values[..., 7])
            rho, omega, g = point_fields(case, r, z, solver.time)
            exact = (rho, r*omega, r*g)
            volume = 2*np.pi*r*mesh.dr*mesh.dz*wr*wz/4
            for index in range(3):
                numerator[index] += np.sum(volume*(measured[index]-exact[index])**2)
                denominator[index] += np.sum(volume*exact[index]**2)
    return {name: float(np.sqrt(num/den)) if den > 0 else None for name, num, den in
            zip(('density', 'rotation_velocity', 'toroidal_field'), numerator, denominator)}
