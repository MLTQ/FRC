"""Independent analytic parent-equation terms at a fixed Fourier/Bessel mode."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from scipy.special import j0,j1


def terms(mode,r,z):
    r,z=np.broadcast_arrays(r,z);ref=mode.ref;q=mode.q
    a=mode.amplitudes(z);az=mode.amplitudes(z,1);azz=mode.amplitudes(z,2)
    f=ref.fields(z);n=f['number_density'];k=2*np.pi/ref.length
    nz=-ref.number_density*ref.density_amplitude*k*np.sin(k*z)
    pez=ref.number_density*E*ref.electron_ev*ref.pressure_amplitude*k*np.cos(k*z+.23)
    p=np.stack((f['electron_pressure'],f['ion_pressure']),axis=-1);pz=np.stack((pez,-pez),axis=-1)
    gamma=np.array([ref.gamma_e,ref.gamma_i]);mass=ref.ion_mass+ref.electron_mass
    alpha=ref.electron_mass/mass;beta=ref.ion_mass/mass;mu=ref.ion_mass*alpha;fractions=np.array([-beta,alpha])
    dp=p*(gamma*a[...,:1]/n[...,None]+a[...,4:6])
    dpz=pz*(gamma*a[...,:1]/n[...,None]+a[...,4:6])
    dpz+=p*(gamma*(az[...,:1]/n[...,None]-a[...,:1]*nz[...,None]/n[...,None]**2)+az[...,4:6])
    j=np.stack((-az[...,7],az[...,6]+q*a[...,8],q*a[...,7]),axis=-1)/MU0
    jz=np.stack((-azz[...,7],azz[...,6]+q*az[...,8],q*az[...,7]),axis=-1)/MU0
    j0r=j0(q*r);j1r=j1(q*r);basis=np.stack((j1r,j1r,j0r),axis=-1)
    electric={};induction={}
    def add(name,coeff,derivative):
        electric[name]=coeff*basis
        induction[name]=np.stack((derivative[...,1]*j1r,-(derivative[...,0]+q*coeff[...,2])*j1r,-q*coeff[...,1]*j0r),axis=-1)
    zero=np.zeros_like(n)
    add('ideal',ref.bias*np.stack((-a[...,2],a[...,1],zero),axis=-1),ref.bias*np.stack((-az[...,2],az[...,1],zero),axis=-1))
    factor=(beta-alpha)*ref.bias/(E*n)
    coeff=factor[...,None]*np.stack((j[...,1],-j[...,0],zero),axis=-1)
    deriv=factor[...,None]*(np.stack((jz[...,1],-jz[...,0],zero),axis=-1)-nz[...,None]/n[...,None]*np.stack((j[...,1],-j[...,0],zero),axis=-1))
    add('hall',coeff,deriv)
    factor=mu*mode.growth/(E**2*n)
    add('inertia',factor[...,None]*j,factor[...,None]*(jz-nz[...,None]/n[...,None]*j))
    radial=-q*np.sum(fractions*dp,axis=-1)/(E*n)
    axial=np.sum(fractions*(dpz/n[...,None]-a[...,:1]*pz/n[...,None]**2),axis=-1)/E
    electric['pressure']=np.stack((radial*j1r,zero,axial*j0r),axis=-1)
    toroidal=q*np.sum(fractions*(pz*a[...,:1]-nz[...,None]*dp),axis=-1)/(E*n*n)
    induction['pressure']=np.stack((zero,toroidal*j1r,zero),axis=-1)
    momentum=dict(pressure=np.stack((q*np.sum(dp,axis=-1)*j1r,zero,-np.sum(dpz,axis=-1)*j0r),axis=-1)/(mass*n[...,None]),
                  lorentz=ref.bias*np.stack((j[...,1]*j1r,-j[...,0]*j1r,zero),axis=-1)/(mass*n[...,None]))
    return dict(electric=electric,induction=induction,momentum=momentum)
