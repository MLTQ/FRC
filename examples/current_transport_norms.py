"""Independent continuous norms and physical power for saved transport states."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .current_transport_sampling import TransportSampler
from .current_transport_reference import exact


def measure(d,mesh,ref,arrays,delta,order=16):
    s=TransportSampler(d,mesh,order);weight=s.weight;u=arrays['state'];r=s.R;z=s.Z
    n,v,vg,nt,_=s.fluid(u,arrays['fluid_rate']);j,jg=s.electric(arrays['current'],True)
    target=exact(ref,r,z,delta);metrics={}
    def norms(value,wanted):
        error=np.sum((value-wanted)**2,axis=-1);norm=np.sum(wanted*wanted,axis=-1)
        def region(rows):return float(np.sqrt(np.sum(weight[rows]*error[rows])/np.sum(weight[rows]*norm[rows])))
        rms=np.sqrt(np.sum(weight*norm)/np.sum(weight))
        return dict(relative_L2=region(slice(None)),first_fluid_element=region(s.r<d.radius/d.elements),last_fluid_element=region(s.r>d.radius*(1-1/d.elements)),maximum_over_rms=float(np.sqrt(error.max())/rms))
    powers={}
    for key in ('advection','stretching','self_drift'):
        value=s.electric(arrays[key+'_electric']);metrics[key]=norms(value,target[key]);powers[key+'_work']=float(np.sum(weight*np.sum(j*value,axis=-1)))
    for key in ('stress_rate','omitted_stress_rate'):
        _,_,_,_,vt=s.fluid(u,arrays[key]);force=d.mass*n[...,None]*vt
        metrics[key]=norms(force,target['force'])
        for k,name in enumerate(('radial','azimuthal','axial')):metrics[key+'_'+name]=norms(force[...,k,None],target['force'][...,k,None])
        powers[key+'_physical_work']=float(np.sum(weight*np.sum(v*force,axis=-1)))
    density=ref.reduced_mass/(2*E**2)*np.sum(j*j,axis=-1)*nt/n**2
    jt=s.electric(arrays['current_rate']);b=s.magnetic(arrays['magnetic']);b[...,2]+=ref.axial_bias;bt=s.magnetic(arrays['magnetic_rate'])
    current=ref.reduced_mass/E**2*(np.sum(j*jt,axis=-1)/n-.5*np.sum(j*j,axis=-1)*nt/n**2)
    magnetic=np.sum(b*bt,axis=-1)/MU0
    powers.update(density_exchange=float(np.sum(weight*density)),current_energy_rate=float(np.sum(weight*current)),magnetic_energy_rate=float(np.sum(weight*magnetic)))
    divj=np.trace(jg,axis1=-2,axis2=-1)
    continuity=nt+n*np.trace(vg,axis1=-2,axis2=-1)+v[...,0]*s.scalar(u[...,0],axis=0)+v[...,2]*s.scalar(u[...,0],axis=2)
    defects=dict(continuity_relative=float(np.sqrt(np.sum(weight*continuity**2)/np.sum(weight*ref.fields(r,z)['nt']**2))),charge_divergence_scaled=float(np.sqrt(ref.radius**2*np.sum(weight*divj**2)/np.sum(weight*np.sum(j*j,axis=-1)))))
    # Actual continuous fluid energy derivative, distinct from grid work identity.
    rate=arrays['full_fluid_rate'];_,_,_,nt,vt=s.fluid(u,rate)
    bulk=d.mass*(.5*nt*np.sum(v*v,axis=-1)+n*np.sum(v*vt,axis=-1));thermal=np.zeros_like(n)
    for k,gamma in enumerate(d.thermal.gammas):
        c=s.scalar(u[...,k+4]);ct=s.scalar(rate[...,k+4]);entropy=c/n;internal=n*E*np.exp(entropy)*(n/d.thermal.n0)**(gamma-1)/(gamma-1)
        thermal+=internal*((gamma-entropy)*nt+ct)/n
    total=bulk+thermal+current+magnetic;scale=float(np.sum(weight*(abs(bulk)+abs(thermal)+abs(current)+abs(magnetic))))
    powers.update(continuous_bulk_rate=float(np.sum(weight*bulk)),continuous_thermal_rate=float(np.sum(weight*thermal)),continuous_total_rate=float(np.sum(weight*total)),continuous_total_backward=float(abs(np.sum(weight*total))/scale),continuous_rate_scale=scale)
    return dict(norms=metrics,powers=powers,defects=defects)
