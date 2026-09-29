"""Independent local pressure/enthalpy energy derivatives and face integrals."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0


def fields(reference, r, z, *, barotropic=False):
    f = reference.fields(r, z)
    if barotropic:
        s, k = (r/reference.radius)**2, 2*np.pi/reference.length
        ln_r_over_r = -4*reference.density_contrast/reference.radius**2*(1-s)*np.cos(k*z)
        ln_z = -reference.density_contrast*k*(1-s)**2*np.sin(k*z)
        p = reference.density*reference.temperature_ev*E*(f['density']/reference.density)**reference.gamma
        t = p/f['density']
        f.update(pressure=p, temperature_j=t, er_over_r=-reference.gamma*t*ln_r_over_r/E,
                 ez=-reference.gamma*t*ln_z/E, g_rate=np.zeros_like(p), entropy_rate=np.zeros_like(p))
    er, ez, b = r*f['er_over_r'], f['ez'], r*f['g']
    enthalpy_factor = -reference.gamma*f['temperature_j']/((reference.gamma-1)*E)
    poynting_r, poynting_z = -ez*b/MU0, er*b/MU0
    enthalpy_r, enthalpy_z = enthalpy_factor*f['jr'], enthalpy_factor*f['jz']
    thermal = (f['density']/reference.density)**(reference.gamma-1)/(reference.gamma-1)*f['entropy_rate']
    return f | dict(energy_rate=thermal+r*r*f['g']*f['g_rate']/MU0,
                    pressure_work=f['jr']*er+f['jz']*ez,
                    radial_energy_flux=poynting_r+enthalpy_r, axial_energy_flux=poynting_z+enthalpy_z)


def cell_average(reference, mesh, key, *, order=12, barotropic=False):
    x, w = np.polynomial.legendre.leggauss(order)
    lo, hi = .5*mesh.rf[:-1]**2, .5*mesh.rf[1:]**2
    result = np.zeros((mesh.nr, mesh.nz))
    for a, wa in zip(x, w):
        r = np.sqrt(lo+hi+a*(hi-lo))[:, None]
        for b, wb in zip(x, w):
            f = fields(reference, r, mesh.z[None]+.5*mesh.dz*b, barotropic=barotropic)
            result += .25*wa*wb*f[key]
    return result


def face_integrals(reference, mesh, *, order=12, barotropic=False):
    x, w = np.polynomial.legendre.leggauss(order)
    fr, fz = np.zeros((mesh.nr+1, mesh.nz)), np.zeros((mesh.nr, mesh.nz+1))
    lo, hi = .5*mesh.rf[:-1]**2, .5*mesh.rf[1:]**2
    for a, wa in zip(x, w):
        f = fields(reference, mesh.rf[:, None], mesh.z[None]+.5*mesh.dz*a, barotropic=barotropic)
        fr += np.pi*mesh.rf[:, None]*mesh.dz*wa*f['radial_energy_flux']
        f = fields(reference, np.sqrt(lo+hi+a*(hi-lo))[:, None], mesh.zf[None], barotropic=barotropic)
        fz += np.pi*mesh.radial_measure[:, None]*wa*f['axial_energy_flux']
    return fr, fz
