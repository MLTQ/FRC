"""Independent physical sphere quadrature and exact exterior dipole energy."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0 as MU0, electron_mass as ME, elementary_charge as E


def integrate(ref,kind,time=0.,order=80):
    nodes,weights=leggauss(order);s=(nodes[:,None]+1)*ref.radius/2;u=nodes[None,:]
    r=s*np.sqrt(1-u*u);z=s*u;w=2*np.pi*s*s*(ref.radius/2)*weights[:,None]*weights[None,:]
    fields=ref.fields(r,z) if kind=='screening' else ref.fields(r,z,time)
    bias=ref.bias if kind=='screening' else 0.
    benergy=np.sum(w*(fields['Br']**2+(fields['Bz']-bias)**2))/(2*MU0)
    boundary=ref.fields(ref.radius,0.)['a'] if kind=='screening' else ref.fields(ref.radius,0.,time)['a']
    benergy+=4*np.pi*ref.radius**3*(boundary-bias/2)**2/(3*MU0)
    density=ref.density(r,z) if kind=='profiled' else ref.electron_density
    kinetic=0. if density is None else np.sum(w*ME*fields['Jtheta']**2/(2*density*E**2))
    heat=0. if kind=='screening' else np.sum(w*ref.heat_density(r,z,time))
    return {'magnetic_energy_j':float(benergy),'current_kinetic_energy_j':float(kinetic),'ohmic_heat_j':float(heat)}
