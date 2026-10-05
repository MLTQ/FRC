"""Original physical-energy second variation about a motionless plasma."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0


def thermal_quadratic(internal_energy, relative_density, specific_entropy, gamma):
    """Half the conservative (n,C=n*s) thermal Hessian contraction."""
    values=np.broadcast_arrays(internal_energy,relative_density,specific_entropy,gamma)
    u,a,b,g=values
    if not all(np.isfinite(v).all() for v in values) or np.any(u<=0) or np.any(g<=1):
        raise ValueError('Finite positive thermal energy and gamma greater than one required')
    return .5*u*((b+(g-1)*a)**2+(g-1)*a*a)


def densities(number_density,pressure,gammas,mass,reduced_mass,physical,current):
    """Energy in a declared linear physical perturbation, without a bias term."""
    n=np.asarray(number_density);p=np.asarray(pressure);v=np.asarray(physical);current=np.asarray(current)
    if (np.any(n<=0) or not all(np.isfinite(x).all() for x in (n,v,current))
        or not np.isfinite([mass,reduced_mass]).all() or min(mass,reduced_mass)<=0):
        raise ValueError('Positive physical background density and masses required')
    thermal=thermal_quadratic(p/(gammas-1),v[..., :1]/n[..., None],v[..., 4:6],gammas)
    return dict(bulk_kinetic=.5*mass*n*np.sum(v[..., 1:4]**2,axis=-1),
        electron_internal=thermal[..., 0],ion_internal=thermal[..., 1],
        current=.5*reduced_mass/E**2*np.sum(current*current,axis=-1)/n,
        magnetic=.5/MU0*np.sum(v[..., 6:9]**2,axis=-1))
