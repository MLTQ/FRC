"""Independent nonlinear convective-current inertia in the original toroidal field."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0


def fields(ref,r,z,reduced_mass,delta):
    """Analytic derivatives, holomorphic in r,z for independent complex-step curl."""
    r,z=np.broadcast_arrays(r,z);s=(r/ref.radius)**2;k=2*np.pi/ref.length
    h=(1-s)**3;hs=-3*(1-s)**2;hss=6*(1-s)
    phase=k*z+.37;sn=np.sin(phase);cs=np.cos(phase)
    n=ref.density*np.exp(ref.density_contrast*(1-s)**2*np.cos(k*z))
    ar=-4*ref.density_contrast*r/ref.radius**2*(1-s)*np.cos(k*z)
    az=-ref.density_contrast*k*(1-s)**2*np.sin(k*z)
    coefficient=ref.toroidal_slope/(MU0*E*n)
    wr=coefficient*r*k*h*sn;wz=2*coefficient*(h+s*hs)*cs
    wr_r=coefficient*k*((h+2*s*hs)*sn-r*h*sn*ar)
    wr_z=coefficient*r*k*h*(k*cs-sn*az)
    wz_r=coefficient*(4*r/ref.radius**2*(2*hs+s*hss)*cs-2*(h+s*hs)*cs*ar)
    wz_z=2*coefficient*(h+s*hs)*(-k*sn-cs*az)
    curl=wr_z-wz_r
    correction=reduced_mass/E*curl
    er=-delta*wz*correction;ez=delta*wr*correction
    advr=wr*wr_r+wz*wr_z;advz=wr*wz_r+wz*wz_z
    return dict(density=n,drift=np.stack((wr,np.zeros_like(wr),wz),axis=-1),
        star_correction=np.stack((np.zeros_like(correction),correction,np.zeros_like(correction)),axis=-1),
        reduced_convective_electric=np.stack((er,np.zeros_like(er),ez),axis=-1),
        full_convective_electric=-delta*reduced_mass/E*np.stack((advr,np.zeros_like(advr),advz),axis=-1),
        kinetic_gradient=reduced_mass*np.stack((wr*wr_r+wz*wz_r,np.zeros_like(wr),wr*wr_z+wz*wz_z),axis=-1))


def g_rate(ref,r,z,reduced_mass,delta,*,epsilon=1e-20,full=False):
    key='full_convective_electric' if full else 'reduced_convective_electric'
    dr=fields(ref,np.asarray(r,dtype=complex)+1j*epsilon,z,reduced_mass,delta)[key].imag/epsilon
    dz=fields(ref,r,np.asarray(z,dtype=complex)+1j*epsilon,reduced_mass,delta)[key].imag/epsilon
    return (dr[...,2]-dz[...,0])/r
