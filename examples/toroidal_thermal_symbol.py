"""Continuum pressure, joint-current and moving-ion toroidal principal matrices."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0, atomic_mass
from .pressure_split_symbol import symbol as pressure_symbol


def coefficients(ref, r, z, angle):
    f = ref.fields(r, z); n, pe, b = f['density'], f['pressure'], r*f['g']
    s, k = (r/ref.radius)**2, 2*np.pi/ref.length
    ar = -4*ref.density_contrast*r*(1-s)*np.cos(k*z)/ref.radius**2
    az = -ref.density_contrast*k*(1-s)**2*np.sin(k*z)
    kr, kz = np.cos(angle), np.sin(angle)
    a, d, c = pressure_symbol(ref, r, z, angle)
    h = b*((ar+2/r)*kz-az*kr)/(MU0*n*E)
    pr, pz = -n*E*r*f['er_over_r'], -n*E*f['ez']
    f_n = (kr*pz-kz*pr)/(E*n*n)
    return np.broadcast_arrays(n, pe, b, a, d, c, h, f_n)


def matrices(ref, r, z, angle, *, ion_pressure=16., ion_mass=39.948*atomic_mass):
    """Return physical first-derivative matrices; time generator is i |k| M."""
    n, pe, b, a, d, c, h, fn = coefficients(ref, r, z, angle)
    shape = n.shape
    pressure = np.zeros(shape+(2,2)); pressure[...,0,1]=a; pressure[...,1,0]=c; pressure[...,1,1]=d
    joint = pressure.copy(); joint[...,0,0]=h
    # Moving massless-electron system: [delta n, delta v_parallel, delta pi, delta pe, delta Bt].
    rho = n*ion_mass; gamma=ref.gamma
    full=np.zeros(shape+(5,5));full[...,0,1]=-n
    full[...,1,2]=-1/rho;full[...,1,3]=-1/rho;full[...,1,4]=-b/(MU0*rho)
    full[...,2,1]=-gamma*ion_pressure
    full[...,3,0]=-gamma*pe/n*d;full[...,3,1]=-gamma*pe;full[...,3,3]=d;full[...,3,4]=c
    full[...,4,0]=-b/n*d+fn;full[...,4,1]=-b;full[...,4,3]=a;full[...,4,4]=h
    cs=np.sqrt(gamma*(pe+ion_pressure)/rho)
    scale=np.stack((n,cs,rho*cs*cs,rho*cs*cs,np.sqrt(MU0*rho)*cs),axis=-1)
    return pressure,joint,full,scale


def spectrum(ref, r, z, angle, **kwargs):
    pressure,joint,full,scale=matrices(ref,r,z,angle,**kwargs)
    # Similarity scaling is a change of units, not a change of physical mass.
    return dict(pressure_growth_per_k=.5*np.sqrt(np.maximum(0,-(pressure[...,1,1]**2+4*pressure[...,0,1]*pressure[...,1,0]))),
                joint_growth_per_k=.5*np.sqrt(np.maximum(0,-((joint[...,0,0]-joint[...,1,1])**2+4*joint[...,0,1]*joint[...,1,0]))),
                moving_eigenvalues=np.linalg.eigvals(full*scale[...,None,:]/scale[..., :,None]))
