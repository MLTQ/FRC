"""Actual q16/q20 species-rate differences on a common physical grid."""
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .current_transport_sampling import TransportSampler
from .trajectory_density_prefix import fingerprint


def calculate(records):
    mode=ContinuumContactMode();ref=mode.ref
    mesh=MHDMesh(ref.radius,ref.length/2,32,64)
    with np.load('outputs/trajectory_tangent_probe/context.npz') as raw:original_base=raw['state']
    alpha=ref.electron_mass/(ref.ion_mass+ref.electron_mass);beta=ref.ion_mass/(ref.ion_mass+ref.electron_mass)
    rows=[];integrals=[]
    for kind in ('interface','mass','space'):
        directory=Path('outputs')/('species_'+kind+'_audit')
        for elements in ([8,16,32] if kind=='space' else [8]):
            d=RegularCoordinateFluid(ref.radius,ref.length,elements,64,ref.ion_mass+ref.electron_mass,
                                     degree=5,gammas=(ref.gamma_e,ref.gamma_i))
            sample=TransportSampler(d,mesh,20);base=np.broadcast_to(original_base[0,0],d.shape);n=sample.scalar(base[...,0])
            for direction in ('saved_error','independent_mode'):
                stem=(f'F{elements}_' if kind=='space' else '')+direction
                with np.load(directory/(stem+'_q16.npz')) as archive:first=dict(archive)
                with np.load(directory/(stem+'_q20.npz')) as archive:last=dict(archive)
                for label,a in [('electron',-beta),('ion',alpha)]:
                    if kind=='interface':
                        keys=[f'{label}_{projection}_{term}' for projection in ('production','independent')
                              for term in ('advection','volume_charge','transfer','radial_axis','radial_wall','radial_interior','axial_seam','axial_interior','sum','direct')]
                        multiplier=-a/E/d.volume.reshape(-1,d.nz);original=last[label+'_production_Cdot']
                    else:
                        keys=[label+'_'+term for term in ('original','lumped','exact_mass')]
                        multiplier=1.;original=last[label+'_original']
                    full=sample.scalar(original.reshape(d.shape[:-1]))/n
                    full2=float(np.sum(sample.weight*full**2))
                    for key in keys:
                        delta=sample.scalar((multiplier*(first[key]-last[key])).reshape(d.shape[:-1]))/n
                        target=sample.scalar((multiplier*last[key]).reshape(d.shape[:-1]))/n
                        diff2=float(np.sum(sample.weight*delta**2));target2=float(np.sum(sample.weight*target**2))
                        integrals.append([diff2,target2,full2])
                        rows.append(dict(kind=kind,elements=elements,direction=direction,component=key,
                            relative_to_component=float(np.sqrt(diff2/target2)) if target2 else None,
                            relative_to_full=float(np.sqrt(diff2/full2)) if full2 else None,
                            max_absolute=float(np.max(abs(delta))),difference_fingerprint=fingerprint(delta)))
    return rows,np.asarray(integrals)
