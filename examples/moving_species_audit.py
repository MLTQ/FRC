"""Independent parent-equation and local-energy evidence for moving species."""
from dataclasses import replace
import numpy as np
from scipy.constants import electron_mass as ME,elementary_charge as E,mu_0 as MU0
from plasma_bench.finite_mass_moments import SpeciesMasses,moments,reduced_ohm
from .moving_species_reference import MovingSpeciesReference,parent_sources,advection,divergence


def audit():
    records=[];default=MovingSpeciesReference()
    references=(('physical',default),('target',replace(default,density=1e20)),
                ('pair',replace(default,ion_mass=ME)),('massless',replace(default,electron_mass=0.)))
    for label,ref in references:
        r,z=np.meshgrid(ref.radius*np.r_[1e-6,1e-4,np.linspace(.02,.98,15)],np.linspace(-.115,.113,19),indexing='ij')
        t=.23/ref.frequency;a=ref.fields(r,z,t);dr,dz,dt=ref.derivatives(r,z,t)
        forces,heat,continuity=parent_sources(ref,r,z,t);mass=SpeciesMasses(ref.ion_mass,ref.electron_mass)
        grad=lambda key:np.stack((dr[key],np.zeros_like(r),dz[key]),axis=-1)
        conv=lambda x,y:advection(a[x],a[y],dr[y],dz[y],r)
        norm=lambda value:float(np.linalg.norm(value))
        relative=lambda error,scale:norm(error)/max(norm(scale),np.finfo(float).tiny)
        ohm=reduced_ohm(mass,a['n'],a['v'],a['w'],a['b'],grad('pe'),grad('pi'),dt['w'],conv('v','w'),conv('w','v'),conv('w','w'))
        source=(mass.beta*forces['e']-mass.alpha*forces['i'])/(E*a['n'][...,None])
        stress_div=mass.reduced*a['n'][...,None]*conv('w','w')
        bulk=mass.total*a['n'][...,None]*(dt['v']+conv('v','v'))-np.cross(a['j'],a['b'])+grad('pe')+grad('pi')+stress_div
        def curl(key):
            return np.stack((-dz[key][...,1],dz[key][...,0]-dr[key][...,2],dr[key][...,1]+a[key][...,1]/r),axis=-1)
        scalar_r=np.sum(dr['v']*a['w']+a['v']*dr['w']-mass.delta*a['w']*dr['w'],axis=-1)
        scalar_z=np.sum(dz['v']*a['w']+a['v']*dz['w']-mass.delta*a['w']*dz['w'],axis=-1)
        star=a['b']+mass.reduced/E*curl('w')
        transformed=np.cross(mass.delta*a['w']-a['v'],star)-mass.reduced/E*np.cross(a['w'],curl('v'))
        transformed-=mass.reduced/E*a['w']*(dt['n']/a['n'])[...,None]
        transformed+=mass.reduced/(a['n'][...,None]*E**2)*dt['j']+(mass.alpha*grad('pi')-mass.beta*grad('pe'))/(E*a['n'][...,None])
        transformed+=mass.reduced/E*np.stack((scalar_r,np.zeros_like(r),scalar_z),axis=-1)
        m=moments(mass,a['n'],a['v'],a['j'],a['b'],a['pe'],a['pi'],a['e'],gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
        magnetic_rate=np.sum(a['b']*dt['b'],axis=-1)/MU0
        energy_rate=magnetic_rate+dt['pe']/(ref.gamma_e-1)+dt['pi']/(ref.gamma_i-1)
        power=heat['i']+heat['e']
        for key,ms in (('i',ref.ion_mass),('e',ref.electron_mass)):
            energy_rate+=.5*ms*dt['n']*np.sum(a['v'+key]**2,axis=-1)+ms*a['n']*np.sum(a['v'+key]*dt['v'+key],axis=-1)
            power+=np.sum(a['v'+key]*forces[key],axis=-1)
        def mixed_flux(rr,zz):
            f=ref.fields(rr,zz,t)
            return mass.reduced*f['n'][...,None]*np.sum(f['v']*f['w'],axis=-1)[...,None]*f['w']
        def difference(fun,axis,h):
            values=[fun(r+i*h if axis==0 else r,z+i*h if axis==1 else z) for i in (-2,-1,1,2)]
            return (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
        full=lambda rr,zz:ref.species_flux(rr,zz,t)
        balances=[];omitted=[]
        for h in (2e-3*ref.radius,1e-3*ref.radius,5e-4*ref.radius):
            div=divergence(full(r,z),difference(full,0,h),difference(full,1,h),r)
            missing=divergence(mixed_flux(r,z),difference(mixed_flux,0,h),difference(mixed_flux,1,h),r)
            balances.append(relative(energy_rate+div-power,power))
            omitted.append(relative(energy_rate+div-missing-power,power))
        minus_curl=-curl('e')
        metrics=dict(ion_continuity=relative(continuity['i'],dt['n']),electron_continuity=relative(continuity['e'],dt['n']),
            bulk_momentum=relative(bulk-forces['i']-forces['e'],forces['i']+forces['e']),
            ohm_parent=relative(ohm+source-a['e'],a['e']),generalized_field_identity=relative(transformed-ohm,ohm),
            faraday=relative(minus_curl-dt['b'],dt['b']),species_energy_flux=relative(m['total_energy_flux']-full(r,z),full(r,z)))
        assert max(metrics.values())<2e-11,(label,metrics)
        assert balances[-1]<1e-9,(label,balances)
        records.append(dict(name=label,pointwise_relative=metrics,local_energy_refinement=balances,
            omitted_mixed_energy_flux_relative=omitted,
            drift_stress_divergence_max_n_m3=float(np.linalg.norm(stress_div,axis=-1).max()),
            density_time_electric_max_v_m=float(np.linalg.norm(mass.reduced/E*a['w']*(dt['n']/a['n'])[...,None],axis=-1).max()),
            vorticity_electric_max_v_m=float(np.linalg.norm(mass.reduced/E*np.cross(a['w'],curl('v')),axis=-1).max())))
    return records
