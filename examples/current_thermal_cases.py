"""Explicit finite-mass current/thermal verification and control protocol."""
from dataclasses import asdict,replace
from scipy.constants import atomic_mass
from .current_beltrami_reference import CurrentBeltrami
from .electron_reference import ElectronReference


def protocol():
    cases=[]
    def add(name,kind,nr,reference,**changes):
        case=dict(name=name,kind=kind,nr=nr,nz=2*nr,reference=asdict(reference),steps=80,
                  model='two_fluid',ion_mass=39.948*atomic_mass,gamma_e=5/3,gamma_i=5/3,
                  density_scale=1e20,quadrature=5,tolerance=1e-10,max_iterations=80,
                  ion_pressure=16.,end_time=1e-6,phase=1.)
        case.update(changes);cases.append(case)
    ref=ElectronReference()
    for nr in (8,16,32,64,128):
        add(f'original_{nr}','original',nr,ref)
    for steps in (20,40,160,320):
        add(f'original_time_{steps}','original',16,ref,steps=steps)
    for nr in (16,32,64,128):
        add(f'massless_original_{nr}','original',nr,ref,model='hall')
    add('original_quadrature','original',16,ref,quadrature=7)
    add('original_scale','original',16,ref,density_scale=1e19)
    add('original_tolerance','original',16,ref,tolerance=1e-11)
    add('original_fine_half_time','original',128,ref,steps=160)
    wave=CurrentBeltrami(number_density=1e18,amplitude=1e-6)
    for helicity in (-1,1):
        for nr in (8,16,32):
            add(f'wave_h{helicity}_{nr}','wave',nr,replace(wave,helicity=helicity))
    for label,reference in (('target',replace(wave,number_density=1e20)),
                             ('strong',replace(wave,number_density=1e20,amplitude=4e-5)),
                             ('inertial_strong',replace(wave,amplitude=1e-5)),
                             ('massless',replace(wave,model='hall'))):
        for nr in (8,16,32):
            add(f'wave_{label}_{nr}','wave',nr,reference,model=reference.model)
    for steps in (10,20,40,320):
        add(f'wave_time_{steps}','wave',8,wave,steps=steps)
    add('wave_reverse','wave',16,replace(wave,axial_mode=-1))
    add('wave_quadrature','wave',16,wave,quadrature=7)
    return cases
