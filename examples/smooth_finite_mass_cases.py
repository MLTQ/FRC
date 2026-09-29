"""Frozen smooth-transport accuracy protocol and initial-state selection."""
from dataclasses import asdict,replace
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from .finite_mass_rate_measurements import protocol as old_rates,case_state as old_rate_state
from .finite_mass_evolution import protocol as old_trajectories,case_state as old_evolution_state
from .finite_mass_contact_reference import FiniteMassContact
from .finite_mass_contact_case import initialize


def rate_protocol():
    cases=[dict(case,kind='rate',entropy_transport='smooth') for case in old_rates()]
    for branch in (0,2):
        original=next(c for c in cases if c['name']==f'wave_{branch}_32')
        cases.append(dict(original,name=f'tiny_{branch}_32',reference=dict(original['reference'],amplitude=.01)))
    for nz in (16,32,64,128):
        for entropy in ('smooth','limited'):
            cases.append(dict(name=f'contact_{entropy}_{nz}',kind='rate',family='contact',nr=8,nz=nz,
                reference=asdict(FiniteMassContact()),entropy_transport=entropy))
    return cases


def evolution_protocol():
    cases=[dict(case,kind='trajectory',entropy_transport='smooth') for case in old_trajectories()]
    ref=FiniteMassContact(bias=0.)
    def add(name,nz,reference,entropy='smooth',magnetic=False,dt_factor=1.):
        speed=abs(reference.velocity)+np.sqrt((reference.gamma_e*reference.electron_ev+reference.gamma_i*reference.ion_ev)*1.602176634e-19/(reference.ion_mass+reference.electron_mass)/(1-reference.density_amplitude))
        dt=1e-9 if magnetic else .18*reference.length/nz/speed
        cases.append(dict(name=name,kind='trajectory',family='contact',nr=4 if not magnetic else 8,nz=nz,
            reference=asdict(reference),entropy_transport=entropy,dt_s=dt*dt_factor,
            duration_s=5e-8 if magnetic else reference.period))
    for nz in (16,32,64):add(f'contact_{nz}',nz,ref)
    add('contact_limited_32',32,ref,'limited')
    add('contact_reverse_32',32,replace(ref,velocity=-ref.velocity))
    add('contact_half_dt_32',32,ref,dt_factor=.5)
    for nz in (16,32,64):add(f'magnetic_contact_{nz}',nz,replace(ref,bias=.005),magnetic=True)
    return cases


def initial(case):
    if case['family']=='contact':
        ref=FiniteMassContact(**case['reference']);mesh=MHDMesh(ref.radius,ref.length/2,case['nr'],case['nz'])
        state,par=initialize(mesh,ref);return ref,mesh,state,par,0.,None
    if case['kind']=='rate':
        ref,mesh,state,par,source,time=old_rate_state(case);return ref,mesh,state,par,time,source
    return old_evolution_state(case)
