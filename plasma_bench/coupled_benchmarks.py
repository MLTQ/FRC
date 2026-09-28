"""Independent smooth solutions for developmental plasma/vacuum verification."""
from math import comb
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_mesh import MHDMesh
from .geometry import Grid
from .coupled import CoupledParameters, CoupledState, VacuumCoupledMHD
from .coupled_hydro import gas_from_primitive
from .electromagnetics import VacuumField


def aligned_meshes(cells, extent=2):
    """Unit fluid cylinder and a larger aligned field cylinder [SI]."""
    return MHDMesh(1., 1., cells, 2*cells), Grid(float(extent), float(extent), extent*cells+1, 2*extent*cells+1)


def compact_flux(r, z, *, amplitude=.01, width=.65, shift=0.):
    """C3 compact flux; amplitude multiplies sqrt(mu0) in SI units."""
    r, z = np.asarray(r), np.asarray(z)-shift
    return amplitude*np.sqrt(MU0)*r*r*np.maximum(1-(r/width)**2, 0)**4*np.maximum(1-(z/width)**2, 0)**4


def compact_force(r, z, *, amplitude=.01, width=.65, bias=None):
    """Continuum J_theta cross B from analytic polynomial derivatives."""
    r, z = np.broadcast_arrays(r, z)
    bias = np.sqrt(MU0) if bias is None else bias
    a, ap, lap = (np.zeros_like(r) for _ in range(3))
    for k in range(5):
        power, coefficient = 2+2*k, comb(4, k)*(-1)**k/width**(2*k)
        a += coefficient*r**power
        ap += coefficient*power*r**(power-1)
        lap += coefficient*power*(power-2)*r**(power-2)
    inside = r < width
    a, ap, lap = (np.where(inside, value, 0.) for value in (a, ap, lap))
    h = np.maximum(1-(z/width)**2, 0)
    b, bp = h**4, -8*z/width**2*h**3
    bpp = -8/width**2*h**3+48*z*z/width**4*h**2
    scale = amplitude*np.sqrt(MU0)
    current = -scale*(lap*b+a*bpp)/(MU0*r)
    br, bz = -scale*a*bp/r, bias+scale*ap*b/r
    return current*bz, -current*br


def cell_average(mesh, function, order=12):
    """Independent tensor Gauss integral with exact cylindrical volume weights."""
    x, w = np.polynomial.legendre.leggauss(order)
    r = mesh.r[:, None, None, None]+.5*mesh.dr*x[None, None, :, None]
    z = mesh.z[None, :, None, None]+.5*mesh.dz*x[None, None, None, :]
    values = function(r, z)
    weight = r*w[None, None, :, None]*w[None, None, None, :]
    def integrate(value):
        return np.sum(value*weight, axis=(-2, -1))*(mesh.dr/4)/mesh.radial_measure[:, None]
    return tuple(integrate(v) for v in values) if isinstance(values, tuple) else integrate(values)


def passive_advection(cells, *, extent=2, speed=.5, amplitude=1e-6, cfl=.2):
    """Weak compact magnetic packet; compare before reflecting-end waves arrive."""
    mesh, grid = aligned_meshes(cells, extent)
    initial = compact_flux(grid.r[:, None], grid.z[None, :], amplitude=amplitude, width=.5)
    solver = VacuumCoupledMHD(mesh, VacuumField(grid, (), bias=0.),
                             CoupledState(gas_from_primitive(mesh, 1., 1., vz=speed), initial),
                             CoupledParameters(cfl=cfl))
    return solver


def gaussian_diffusion(r, z, time, *, diffusivity=.05, width=.25, amplitude=1e-6):
    """Unbounded analytic solution of psi_t=D Delta*psi with regular axis."""
    initial_width = width**2
    spread = initial_width+4*diffusivity*time
    return amplitude*np.sqrt(MU0)*r*r*(initial_width/spread)**2.5*np.exp(-(r*r+z*z)/spread)


def passive_diffusion(cells, *, diffusivity=.05, amplitude=1e-6, cfl=.2):
    """Weak Gaussian; harmonic exterior differs from the unbounded tail."""
    from .coupling_mesh import CouplingMesh
    mesh, grid = aligned_meshes(cells)
    initial = gaussian_diffusion(grid.r[:, None], grid.z[None, :], 0., diffusivity=diffusivity, amplitude=amplitude)
    initial = CouplingMesh(mesh, grid).harmonic_extension(initial)
    return VacuumCoupledMHD(mesh, VacuumField(grid, (), bias=0.),
                           CoupledState(gas_from_primitive(mesh, 1., 1.), initial),
                           CoupledParameters(resistivity=diffusivity*MU0, cfl=cfl))


def coil_pulse(cells, *, extent=2, cfl=.2):
    """Two prescribed loops drive initially resting chamber gas; no FRC claim."""
    from .coils import Coil
    from .waveforms import Unipolar
    from .electromagnetics import loop_field
    mesh, grid = aligned_meshes(cells, extent)
    radius, positions = 1.23, (-.37, .37)
    unit_center = sum(float(loop_field(0., 0., radius, z)[2]) for z in positions)
    wave = Unipolar(.3*np.sqrt(MU0)/unit_center, rise=.04, hold=.02, fall=.04)
    coils = tuple(Coil(radius, z, wave, name=f"verification_{i}") for i, z in enumerate(positions))
    return VacuumCoupledMHD(mesh, VacuumField(grid, coils, bias=np.sqrt(MU0)),
                           CoupledState(gas_from_primitive(mesh, 1., 1.), np.zeros(grid.shape)),
                           CoupledParameters(resistivity=.1*MU0, cfl=cfl, max_timestep=.002))
