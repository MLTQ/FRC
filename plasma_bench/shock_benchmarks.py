"""Embed an exact local cylindrical magnetic shock in the actual 2-D solvers."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .shock_reference import CylindricalShock
from .mhd_mesh import MHDMesh
from .geometry import Grid
from .coupling_mesh import CouplingMesh
from .coupled import VacuumCoupledMHD, CoupledParameters, CoupledState
from .electromagnetics import VacuumField
from .mhd import AxisymmetricMHD, MHDParameters
from .mhd_state import MHDState


def regularized_profile(reference, radius, *, inner=.38, outer=.85, wall=1., wall_taper=True):
    """The reference is exact only in the declared annulus; regular axis elsewhere."""
    radius = np.asarray(radius)
    rho, velocity, pressure, field = reference.profile(np.clip(radius, inner, outer))
    velocity = np.where(radius < inner, velocity*radius/inner, velocity)
    if wall_taper:
        coordinate = np.clip((radius-outer)/(wall-outer), 0., 1.)
        velocity *= .5*(1+np.cos(np.pi*coordinate))
    return rho, velocity, pressure, field


def initial_fields(mesh, reference=None, *, wall_taper=True, inner=.38):
    reference = reference or CylindricalShock()
    if not np.isclose(reference.radius/mesh.dr, round(reference.radius/mesh.dr), atol=1e-12):
        raise ValueError("Initial shock must align with a radial cell face")
    x, weight = np.polynomial.legendre.leggauss(12)
    r = mesh.r[:, None]+.5*mesh.dr*x[None, :]
    rho, velocity, pressure, field = regularized_profile(reference, r, wall=mesh.radius, wall_taper=wall_taper, inner=inner)
    average = lambda value: np.sum(value*r*weight, axis=1)*(.5*mesh.dr)/mesh.radial_measure
    gas = np.zeros((mesh.nr, mesh.nz, 4))
    gas[..., 0] = average(rho)[:, None]
    gas[..., 1] = average(rho*velocity)[:, None]
    gas[..., 3] = average(pressure/(reference.gamma-1)+.5*rho*velocity**2)[:, None]
    flux = np.concatenate(([0.], np.cumsum(average(field)*mesh.radial_measure)))
    return gas, np.broadcast_to(flux[:, None], (mesh.nr+1, mesh.nz+1)).copy()


def standing_shock(cells, *, coupled=True, cfl=.2, midpoint_method="picard", shock_diffusion=0., wall_taper=True,
                   inner=.38, half_length=1):
    reference = CylindricalShock()
    if half_length not in (1, 2):
        raise ValueError("Shock fixture supports chamber half-lengths 1 or 2 m")
    mesh = MHDMesh(1., float(half_length), cells, 2*half_length*cells)
    extent = half_length+1
    grid = Grid(float(extent), float(extent), extent*cells+1, 2*extent*cells+1)
    gas, total = initial_fields(mesh, reference, wall_taper=wall_taper, inner=inner)
    if not coupled:
        conserved = gas.copy()
        br, bz = mesh.cell_fields(total)
        conserved[..., 3] += (br*br+bz*bz)/(2*MU0)
        return AxisymmetricMHD(mesh, MHDState(conserved, total), MHDParameters(cfl=cfl)), reference
    sources = VacuumField(grid, (), bias=reference.upstream_field)
    mapping = CouplingMesh(mesh, grid, boundary="multipole")
    induced = np.zeros(grid.shape)
    induced[:mesh.nr+1, mapping.z_offset:mapping.z_offset+mesh.nz+1] = total-mapping.fluid_flux(sources.at(0).psi)
    induced = mapping.harmonic_extension(induced)
    solver = VacuumCoupledMHD(mesh, sources, CoupledState(gas, induced),
                             CoupledParameters(cfl=cfl, vacuum_boundary="multipole", midpoint_method=midpoint_method,
                                               shock_diffusion=shock_diffusion))
    return solver, reference
