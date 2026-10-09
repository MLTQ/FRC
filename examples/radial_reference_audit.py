"""Independent physical shape and true-curl controls for radial references."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0
from .archived_contact_mode import ArchivedContactMode
from .matched_mode_rates import fingerprint
from .implicit_plasma_records import provenance,freeze,sha


def curl_controls(mode):
    r=np.array([.009,.017,.029]);z=np.array([-.083,-.014,.071]);records=[];raw={}
    def b(rr,zz):return mode.fields(rr,zz)['physical'][...,6:9]
    target=mode.fields(r,z)['current'];raw['sampled_current']=target
    for i,(hr,hz) in enumerate(((2e-6,2e-5),(1e-6,1e-5),(.5e-6,.5e-5))):
        dr=(b(r-2*hr,z)-8*b(r-hr,z)+8*b(r+hr,z)-b(r+2*hr,z))/(12*hr)
        dz=(b(r,z-2*hz)-8*b(r,z-hz)+8*b(r,z+hz)-b(r,z+2*hz))/(12*hz)
        field=b(r,z);curl=np.stack((-dz[:,1],dz[:,0]-dr[:,2],dr[:,1]+field[:,1]/r),axis=-1)/MU0
        records.append(dict(radial_difference_step=hr,axial_difference_step=hz,
            component_relative_errors=(np.linalg.norm(curl-target,axis=0)/np.linalg.norm(target,axis=0)).tolist()))
        raw[f'finite_difference_curl_{i}']=curl
    return records,raw


def run(directory,destination):
    directory=Path(directory);destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    parent=json.loads((directory/'summary.json').read_text())
    if not parent['terminal']:raise ValueError('Complete reference refinement required')
    modes=[ArchivedContactMode(directory,i) for i in range(len(parent['continuum']))]
    identity=provenance(['examples.radial_reference_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    cases=[]
    for mode in modes:
        curl,raw=curl_controls(mode);quadrature=[]
        for order in (12,16):
            q,w=np.polynomial.legendre.leggauss(order);radius=mode.ref.radius;length=mode.ref.length
            r=(radius/16*(np.arange(16)[:,None]+.5+.5*q)).ravel()
            z=(-length/2+length/32*(np.arange(32)[:,None]+.5+.5*q)).ravel()
            wr=(np.pi*radius/16*r.reshape(16,order)*w).ravel()
            wz=np.broadcast_to(length/64*w,(32,order)).ravel();weight=wr[:,None]*wz[None,:]
            value=mode.fields(r[:,None],z[None,:]);target=modes[-1].fields(r[:,None],z[None,:])
            # Eigenvector sign and overall amplitude are arbitrary: measure the
            # density-based gauge explicitly, never alter an evolution state.
            a=value['physical'][...,0];b=target['physical'][...,0]
            gain=float(np.sum(weight*a*b)/np.sum(weight*a*a));errors={}
            for name in ('physical','current','electric'):
                error=gain*value[name]-target[name]
                errors[name]=(np.sqrt(np.sum(weight[...,None]*error**2,axis=(0,1))/
                    np.sum(weight[...,None]*target[name]**2,axis=(0,1)))).tolist()
            quadrature.append(dict(order=order,density_alignment_gain=gain,component_relative_errors=errors,
                physical_fingerprints={name:fingerprint(v) for name,v in value.items()}))
        raw['physical_nodes']=mode.physical_nodes;raw['fourier_coefficients']=mode.coefficients
        path=destination/f'case_{len(cases):02}.npz';np.savez_compressed(path,**raw)
        cases.append(dict(points=mode.points,radial_mode=mode.radial_mode,continuum_growth=mode.growth,input_sha256=mode.input_sha256,
            curl_controls=curl,quadrature=quadrature,nyquist_coefficient_fraction=(abs(mode.coefficients[:,mode.points//2])/
                np.max(abs(mode.coefficients),axis=1)).tolist(),arrays_file=path.name,arrays_sha256=sha(path)))
    if provenance(['examples.radial_reference_audit'])!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_summary_sha256=sha(directory/'summary.json'),
        target_points=modes[-1].points,cases=cases,
        interpretation='Finite-difference curl includes the real Nyquist cosine derivative; the original collocation derivative sets that component to zero. Retain the mismatch, refine the reference, and never suppress the coefficient.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('complete independent reference physical audit',destination,flush=True)
    return result
