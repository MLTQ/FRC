"""Matched continuum states and initialization for physical log-thermal checks."""
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.compatible_log_thermal import CompatibleLogThermal
from plasma_bench.entropy_scalar import EntropyScalar
from .variable_current_verification import initialize as current_initialize


def target(reference,r,z,gamma,kind='smooth',density_scale=1e20):
    f=reference.fields(r,z);n=f['n']
    if kind=='smooth':
        thermal=reference.thermal(r,z,gamma=gamma,density_scale=density_scale)
        s=np.log(thermal['entropy']/n/E)
        st=thermal['entropy_rate']/thermal['entropy']
        p=thermal['pressure'];electric=thermal['electric']
    elif kind=='barotropic':
        s=np.full_like(n,2.);st=np.zeros_like(n)
        p=n*E*(n/density_scale)**(gamma-1)*np.exp(s)
        electric=-gamma*(p/(E*n*n))[...,None]*f['gradn']
    elif kind=='isothermal':
        temperature=2.;s=np.log(temperature)-(gamma-1)*np.log(n/density_scale)
        st=-(gamma-1)*np.sum(f['current']*f['gradn'],axis=-1)/(E*n*n)
        p=n*E*temperature;electric=-temperature*f['gradn']/n[...,None]
    else:raise ValueError('Unknown thermal reference kind')
    proxy=n*E*np.exp(s);internal=p/(gamma-1)
    return dict(entropy=s,entropy_rate=st,pressure=p,electric=electric,proxy=proxy,
        proxy_rate=proxy*st,internal=internal,internal_entropy_rate=internal*st,
        internal_density_rate=gamma*internal*f['nt']/n,temperature=p/(E*n))


def initialize(spec):
    physical=dict(spec)
    if spec.get('kind','smooth')!='smooth':physical['speed']=0.
    reference,model,field,nt,radial,axial=current_initialize(physical)
    if spec.get('direction',1)==-1:
        reference=replace(reference,amplitude=-reference.amplitude,toroidal=-reference.toroidal)
        field=HarmonicMagnetic(model.mesh,-field.psi_increment,-field.g,field.axial_bias)
    scalar=EntropyScalar(model.density_map)
    species={label:CompatibleLogThermal(model,gamma=gamma,scalar=scalar)
             for label,gamma in (('electron',5/3),('ion',1.4))}
    r,z=np.meshgrid(model.mesh.rf,model.mesh.zf[:-1],indexing='ij')
    states={label:target(reference,r,z,s.gamma,spec.get('kind','smooth'))['entropy'] for label,s in species.items()}
    return reference,model,field,nt,radial,axial,species,states


def protocol():
    specs=[dict(nr=n,degree=p) for p in (3,5) for n in (8,16,32,64)]
    specs.extend(dict(nr=32,quadrature_order=q) for q in (8,12))
    specs.extend(dict(nr=32,number_density=n) for n in (1e16,1e20))
    specs.append(dict(nr=32,direction=-1))
    specs.append(dict(nr=16,kind='barotropic'))
    specs.extend(dict(nr=n,kind='isothermal') for n in (8,16,32,64))
    return specs
