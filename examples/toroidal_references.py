"""Independent exact equilibria/transport and asymptotic toroidal wave fixtures."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.special import j0, j1, jn_zeros
from plasma_bench.toroidal_state import ToroidalState, from_primitive


def r2_average(mesh):
    return ((mesh.rf[1:]**4-mesh.rf[:-1]**4)/(2*(mesh.rf[1:]**2-mesh.rf[:-1]**2)))[:, None]


def equilibrium(mesh, *, omega=.4, current_slope=.3, pressure=1., balanced=True, gamma=5/3):
    """rho=1, vtheta=omega r, Btheta=sqrt(mu0) current_slope r."""
    r, _ = mesh.corners
    p = pressure+(.5*omega**2-current_slope**2)*r2_average(mesh) if balanced else pressure
    return from_primitive(mesh, .5*np.sqrt(MU0)*r*r, 1., p, omega=omega,
                          toroidal_over_radius=np.sqrt(MU0)*current_slope, gamma=gamma)


def rotating_contact(mesh, time=0., *, contrast=.5, omega=.4, current_slope=.3, velocity=.8, gamma=5/3):
    """Exact finite-amplitude translating rotating entropy wave, Bz=Br=0."""
    x, weight = np.polynomial.legendre.leggauss(12)
    z = mesh.z[:, None]+mesh.dz*x/2
    density = 1+contrast*np.cos(np.pi*(z-velocity*time)/mesh.half_length)
    rho = .5*np.sum(density*weight, axis=-1)[None, :]
    angular_coefficient = omega*.5*np.sum(np.sqrt(density)*weight, axis=-1)[None, :]
    r2 = r2_average(mesh)
    p = 1+(.5*omega**2-current_slope**2)*r2
    energy = p/(gamma-1)+.5*rho*velocity**2+.5*r2*(omega**2+current_slope**2)
    shape = (mesh.nr, mesh.nz)
    fields = (rho, 0., rho*velocity, energy, r2*angular_coefficient, np.sqrt(MU0)*current_slope)
    u = np.stack([np.broadcast_to(v, shape) for v in fields], axis=-1).copy()
    return ToroidalState(u, np.zeros((mesh.nr+1, mesh.nz+1)))


def torsional_wave(mesh, time=0., *, amplitude=1e-4, axial_mode=1, gamma=5/3):
    """Linear torsional Alfven wave at vA=rho=p=1, errors O(amplitude²)."""
    k = np.pi*axial_mode/mesh.half_length
    mode = np.sin(k*(mesh.z-time))[None, :]*np.sinc(k*mesh.dz/(2*np.pi))
    r, _ = mesh.corners
    return from_primitive(mesh, .5*np.sqrt(MU0)*r*r, 1., 1., omega=amplitude*mode,
                          toroidal_over_radius=-np.sqrt(MU0)*amplitude*mode, gamma=gamma)


def toroidal_decay(mesh, time=0., *, diffusivity=.1, amplitude=1e-7, gamma=5/3):
    """Btheta=C J1(qr) exp(-Dq²t), qR first J0 zero, Ez(R)=0."""
    k = jn_zeros(0, 1)[0]/mesh.radius
    a, b = mesh.rf[:-1], mesh.rf[1:]
    mean_g = amplitude*(j0(k*a)-j0(k*b))/(k*mesh.radial_measure)*np.exp(-diffusivity*k*k*time)
    state = from_primitive(mesh, np.zeros((mesh.nr+1, mesh.nz+1)), 1., 1., toroidal_over_radius=mean_g[:, None], gamma=gamma)
    return state, k


def decay_field(radius, time, *, mesh, diffusivity=.1, amplitude=1e-7):
    k = jn_zeros(0, 1)[0]/mesh.radius
    return amplitude*j1(k*radius)*np.exp(-diffusivity*k*k*time)


def point_fields(case, r, z, time):
    """Continuum rho, Omega and g at arbitrary physical quadrature points."""
    r, z = np.broadcast_arrays(r, z)
    rho = np.ones_like(r)
    omega = np.full_like(r, case.omega)
    g = np.full_like(r, np.sqrt(MU0)*case.current_slope)
    if case.kind == 'contact':
        rho = 1+case.contrast*np.cos(np.pi*(z-case.velocity*time)/case.half_length)
        omega = case.omega/np.sqrt(rho)
    elif case.kind == 'wave':
        k = np.pi*case.axial_mode/case.half_length
        omega = case.amplitude*np.sin(k*(z-time))
        g = -np.sqrt(MU0)*omega
    elif case.kind == 'decay':
        k = jn_zeros(0, 1)[0]/case.radius
        ratio = np.divide(j1(k*r), r, out=np.full_like(r, k/2), where=r != 0)
        g = case.amplitude*ratio*np.exp(-case.diffusivity*k*k*time)
        omega.fill(0.)
    elif case.kind != 'equilibrium':
        raise ValueError('Unknown continuum reference')
    return rho, omega, g
