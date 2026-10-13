"""Independent continuum fluid terms and actual-reconstruction interior terms."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.special import j0, j1


FLUID_NAMES=('density','radial_velocity','azimuthal_velocity','axial_velocity','electron_entropy','ion_entropy')


def magnetic_responses(mode, r, z, normalization):
    ref=mode.ref;count=mode.points;nodes=np.arange(count)*ref.length/count-ref.length/2
    a=mode.physical_nodes/normalization;az=(mode.derivative@a.T).T
    f=ref.fields(nodes);n=f['number_density'];p=np.array([f['electron_pressure'],f['ion_pressure']])
    k=2*np.pi/ref.length;npz=-ref.number_density*ref.density_amplitude*k*np.sin(k*nodes)
    pez=ref.number_density*E*ref.electron_ev*ref.pressure_amplitude*k*np.cos(k*nodes+.23)
    pz=np.array([pez,-pez]);gamma=np.array([ref.gamma_e,ref.gamma_i])[:,None]
    dp=p*(gamma*a[0]/n+a[4:6])
    dpz=pz*(gamma*a[0]/n+a[4:6])+p*(gamma*(az[0]/n-a[0]*npz/n**2)+az[4:6])
    mass=ref.ion_mass+ref.electron_mass;alpha=ref.electron_mass/mass;beta=ref.ion_mass/mass
    mu=ref.ion_mass*alpha;current=(mode.curl@a[6:].ravel()).reshape(3,count)/MU0
    coefficient=mu/(E**2*n)
    inertia=np.eye(3*count)+mode.curl@np.diag(np.tile(coefficient,3))@mode.curl/MU0
    electric={};forcing={}
    ideal=np.vstack((-ref.bias*a[2],ref.bias*a[1],np.zeros(count)))
    hall=np.vstack(((beta-alpha)*ref.bias*current[1]/(E*n),-(beta-alpha)*ref.bias*current[0]/(E*n),np.zeros(count)))
    for name,value in [('ideal_lorentz',ideal),('hall',hall)]:
        electric[name]=value;forcing[name]=-mode.curl@value.ravel()
    for i,(name,fraction) in enumerate(zip(('electron_species','ion_species'),(-beta,alpha))):
        electric[name]=fraction/E*np.vstack((-mode.q*dp[i]/n,np.zeros(count),dpz[i]/n-a[0]*pz[i]/n**2))
        rhs=np.zeros((3,count));rhs[1]=fraction*mode.q*(pz[i]*a[0]-npz*dp[i])/(E*n**2)
        forcing[name]=rhs.ravel()
    r,z=np.broadcast_arrays(r,z);unique,inverse=np.unique(z,return_inverse=True)
    phase=np.exp(1j*(unique[:,None]+ref.length/2)*mode.frequency[None,:])
    phase[:,count//2]=np.cos(np.pi*count*(unique+ref.length/2)/ref.length)
    basis=np.stack((j1(mode.q*r),j1(mode.q*r),j0(mode.q*r)),axis=-1)
    def expand(value):
        values=(phase@(np.fft.fft(value,axis=1)/count).T).real
        return values[inverse].reshape(z.shape+(3,))*basis
    result={}
    for name,rhs in forcing.items():
        bt=np.linalg.solve(inertia,rhs).reshape(3,count)
        jt=(mode.curl@bt.ravel()).reshape(3,count)/MU0
        et=electric[name]+coefficient*jt
        result[name]=(expand(bt),expand(jt),expand(et))
    return result


def named(rates, n, pressure, gammas):
    result={key:rates[...,i] for i,key in enumerate(FLUID_NAMES)}
    pt=pressure*(gammas*rates[...,:1]/n[...,None]+rates[...,4:6])
    result.update(electron_pressure=pt[...,0],ion_pressure=pt[...,1])
    return result


def analytic(mode, r, z, normalization):
    r,z=np.broadcast_arrays(r,z);ref=mode.ref;f=ref.fields(z)
    a=mode.amplitudes(z)/normalization;az=mode.amplitudes(z,1)/normalization
    n=f['number_density'];p=np.stack((f['electron_pressure'],f['ion_pressure']),axis=-1)
    k=2*np.pi/ref.length;npz=-ref.number_density*ref.density_amplitude*k*np.sin(k*z)
    pez=ref.number_density*E*ref.electron_ev*ref.pressure_amplitude*k*np.cos(k*z+.23)
    pz=np.stack((pez,-pez),axis=-1);gamma=np.array([ref.gamma_e,ref.gamma_i])
    sgrad=pz/p-gamma*npz[...,None]/n[...,None]
    dp=p*(gamma*a[...,:1]/n[...,None]+a[...,4:6])
    dpz=pz*(gamma*a[...,:1]/n[...,None]+a[...,4:6])
    dpz+=p*(gamma*(az[...,:1]/n[...,None]-a[...,:1]*npz[...,None]/n[...,None]**2)+az[...,4:6])
    scalar=j0(mode.q*r);transverse=j1(mode.q*r)
    mass=ref.ion_mass+ref.electron_mass;alpha=ref.electron_mass/mass;beta=ref.ion_mass/mass
    current=mode.fields(r,z)['current']/normalization
    bulk=np.zeros(r.shape+(6,));lorentz=np.zeros_like(bulk)
    bulk[...,0]=(-mode.q*n*a[...,1]-npz*a[...,3]-n*az[...,3])*scalar
    bulk[...,1]=mode.q*dp.sum(-1)*transverse/(mass*n)
    bulk[...,3]=-dpz.sum(-1)*scalar/(mass*n)
    bulk[...,4:6]=-sgrad*(a[...,3]*scalar)[...,None]
    lorentz[...,1]=ref.bias*current[...,1]/(mass*n)
    lorentz[...,2]=-ref.bias*current[...,0]/(mass*n)
    terms=dict(bulk=bulk,ideal_lorentz=lorentz)
    for species,(name,fraction) in enumerate(zip(('electron_species','ion_species'),(-beta,alpha))):
        value=np.zeros_like(bulk);value[...,4+species]=-fraction*sgrad[...,species]*current[...,2]/(E*n)
        terms[name]=value
    result={name:named(value,n,p,gamma) for name,value in terms.items()}
    responses=magnetic_responses(mode,r,z,normalization)
    for name,(bt,jt,et) in responses.items():
        if name not in result:result[name]=named(np.zeros_like(bulk),n,p,gamma)
        for prefix,values in [('B',bt),('current',jt),('electric',et)]:
            for k,axis in enumerate(('radial','azimuthal','axial')):
                key=axis+'_B' if prefix=='B' else prefix+'_'+axis
                result[name][key]=values[...,k]
    return result


def interior(sample, state, du, current, bias, ion_mass, electron_mass):
    if np.any(state[...,1:4]):raise ValueError('Stationary zero-velocity base required')
    d=sample.d;n=sample.scalar(state[...,0]);dn=sample.scalar(du[...,0])
    s=np.stack([sample.scalar(state[...,k])/n for k in (4,5)],axis=-1)
    ds=np.stack([(sample.scalar(du[...,k])-s[...,k-4]*dn)/n for k in (4,5)],axis=-1)
    gamma=d.thermal.gammas;p=E*n[...,None]*np.exp(s)*(n[...,None]/d.thermal.n0)**(gamma-1)
    gradient_n=np.zeros(n.shape+(3,));gradient_dn=np.zeros_like(gradient_n)
    gradient_s=np.zeros(s.shape+(3,));gradient_ds=np.zeros_like(gradient_s)
    for axis in (0,2):
        gn=sample.scalar(state[...,0],axis=axis);gdn=sample.scalar(du[...,0],axis=axis)
        gradient_n[...,axis]=gn;gradient_dn[...,axis]=gdn
        for k in (4,5):
            gs=(sample.scalar(state[...,k],axis=axis)-s[...,k-4]*gn)/n
            gradient_s[...,k-4,axis]=gs
            gradient_ds[...,k-4,axis]=(sample.scalar(du[...,k],axis=axis)-gs*dn-s[...,k-4]*gdn-ds[...,k-4]*gn)/n
    gp=p[...,None]*(gamma[None,None,:,None]*gradient_n[...,None,:]/n[...,None,None]+gradient_s)
    gdp=gp*(gamma*dn[...,None]/n[...,None]+ds)[...,None]
    gdp+=p[...,None]*(gamma[None,None,:,None]*(gradient_dn/n[...,None]-dn[...,None]*gradient_n/n[...,None]**2)[...,None,:]+gradient_ds)
    velocity=np.stack((sample.R*sample.scalar(du[...,1]/state[...,0],True),
                       sample.R*sample.scalar(du[...,2]/state[...,0],True),sample.scalar(du[...,3])/n),axis=-1)
    radial=sample.R*sample.scalar(du[...,1]/state[...,0],True,0)+velocity[...,0]/sample.R
    axial=(sample.scalar(du[...,3],axis=2)-velocity[...,2]*gradient_n[...,2])/n
    divv=radial+velocity[...,0]/sample.R+axial
    j,gradj=sample.electric(current,gradient=True);divj=np.trace(gradj,axis1=-2,axis2=-1)
    mass=ion_mass+electron_mass;alpha=electron_mass/mass;beta=ion_mass/mass
    bulk=np.zeros(n.shape+(6,));lorentz=np.zeros_like(bulk)
    bulk[...,0]=-n*divv-np.sum(velocity*gradient_n,axis=-1)
    bulk[...,1:4]=-gdp.sum(-2)/(mass*n[...,None])+dn[...,None]*gp.sum(-2)/(mass*n[...,None]**2)
    bulk[...,4:6]=-np.sum(gradient_s*velocity[...,None,:],axis=-1)
    lorentz[...,1:4]=np.cross(j,[0.,0.,bias])/(mass*n[...,None])
    terms=dict(bulk=bulk,ideal_lorentz=lorentz);charge={}
    for species,(name,fraction) in enumerate(zip(('electron_species','ion_species'),(-beta,alpha))):
        value=np.zeros_like(bulk);value[...,4+species]=-fraction*np.sum(j*gradient_s[...,species,:],axis=-1)/(E*n)
        terms[name]=value
        extra=np.zeros_like(bulk);extra[...,4+species]=-fraction*s[...,species]*divj/(E*n)
        charge[name]=dict(conservative=named(value+extra,n,p,gamma),skew_volume=named(value+.5*extra,n,p,gamma))
    integral=lambda x:float(np.sum(sample.weight*x))
    j2=integral(np.sum(j*j,axis=-1));div2=integral(divj*divj)
    diagnostics=dict(current_squared_integral=j2,current_divergence_squared_integral=div2,
        radius_times_divergence_over_current=float(d.radius*np.sqrt(div2/j2)) if j2 else None,
        total_base_pressure_gradient_RMS=float(np.sqrt(integral(np.sum(gp.sum(-2)**2,axis=-1))/sample.weight.sum())))
    return {name:named(value,n,p,gamma) for name,value in terms.items()},charge,diagnostics
