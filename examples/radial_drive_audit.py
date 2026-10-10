"""Independent forcing consistency even when full induction is rejected."""
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.strong_wall_plasma import StrongWallFluid
from plasma_bench.mhd_mesh import MHDMesh
from .archived_contact_mode import ArchivedContactMode
from .radial_parent_terms import terms
from .current_transport_sampling import TransportSampler
from .radial_term_sampling import radial_slice,velocity_rate,projected_faraday,wall_projection
from .implicit_plasma_records import provenance,freeze,sha


MOMENTUM=('bulk','lorentz','drift','species','wall_reaction','total')
DRIVE=('ideal','hall','density','transport','pressure_full','pressure_reduced','pressure_potential')


def measure(record,raw,order):
    mode=ArchivedContactMode();ref=mode.ref;sp=record['spec'];epsilon=record['epsilon']
    d=StrongWallFluid(ref.radius,ref.length,sp['elements'],sp['points'],ref.ion_mass+ref.electron_mass,gammas=(ref.gamma_e,ref.gamma_i))
    mesh=MHDMesh(ref.radius,ref.length/2,sp['field_nr'],sp['field_nz']);sample=TransportSampler(d,mesh,order);p=wall_projection(mesh,sp.get('quadrature',12))
    center=lambda k:(raw['plus_'+k]-raw['minus_'+k])/(2*epsilon)
    magnetic={k:projected_faraday(mesh,center('electric_'+k),p) for k in DRIVE}
    magnetic['driving_sum']=sum(magnetic[k] for k in ('ideal','hall','density','transport','pressure_reduced'))
    arrays={};nets={}
    for start in range(0,mesh.nr,8):
        stop=min(start+8,mesh.nr);s=radial_slice(sample,start*order,stop*order);count=stop-start
        analytic=terms(mode,s.R,s.Z);fields=mode.fields(s.R,s.Z)
        integral=lambda a:np.sum((s.weight[...,None]*a).reshape(count,order,len(s.z),3),axis=(1,2))
        for group,target in [('momentum',mode.growth*fields['physical'][...,1:4]),('driving',mode.growth*fields['physical'][...,6:9])]:
            if group not in nets:nets[group]=np.zeros((mesh.nr,3))
            nets[group][start:stop]=integral(target**2)
        def add(group,name,value,target):
            key=group+'_'+name
            if key not in arrays:arrays[key]=np.zeros((mesh.nr,3,3))
            arrays[key][start:stop]=np.stack((integral(value**2),integral(target**2),integral((value-target)**2)),axis=1)
        for name in MOMENTUM:
            value=(velocity_rate(s,raw['plus_state'],raw['plus_fluid_'+name])-velocity_rate(s,raw['minus_state'],raw['minus_fluid_'+name]))/(2*epsilon)
            target=sum(analytic['momentum'].values()) if name=='total' else analytic['momentum'].get('pressure' if name=='bulk' else name,np.zeros_like(value))
            add('momentum',name,value,target)
        for name,b in magnetic.items():
            value=s.magnetic(b)
            if name=='driving_sum':target=sum(analytic['induction'][k] for k in ('ideal','hall','pressure'))
            else:target=analytic['induction'].get('pressure' if name in ('pressure_full','pressure_reduced') else name,np.zeros_like(value))
            add('driving',name,value,target)
    metrics={}
    for key,values in arrays.items():
        group=key.split('_')[0];regions={}
        for name,region in [('volume',slice(None)),('axis',slice(0,mesh.nr//d.elements)),('wall',slice(mesh.nr-mesh.nr//d.elements,None))]:
            a,t,e=values[region].sum(axis=0);n=nets[group][region].sum(axis=0)
            regions[name]=dict(actual_L2=np.sqrt(a).tolist(),target_L2=np.sqrt(t).tolist(),error_L2=np.sqrt(e).tolist(),
                error_over_full_reference_rate_L2=np.sqrt(e/np.maximum(n,np.finfo(float).tiny)).tolist(),
                relative_to_own_term=[float(np.sqrt(x/y)) if y>0 else None for x,y in zip(e,t)])
        metrics[key]=regions
    arrays.update({'net_'+k:v for k,v in nets.items()});arrays['radius']=mesh.r
    return dict(order=order,metrics=metrics),arrays


def run(destination,original='outputs/radial16_drive_capture'):
    destination=Path(destination);root=Path(original)
    if destination.exists():raise ValueError('Refusing overwrite')
    source=json.loads((root/'summary.json').read_text())
    if not source['terminal']:raise ValueError('Terminal input required')
    identity=provenance(['examples.radial_drive_audit']);destination.mkdir(parents=True);freeze(destination,identity);records=[]
    for i,record in enumerate(source['cases']):
        start=time.monotonic();path=root/record['arrays_file']
        if sha(path)!=record['arrays_sha256']:raise ValueError('Input changed')
        raw=dict(np.load(path));quadrature=[];arrays={}
        for order in (16,20):
            result,data=measure(record,raw,order);quadrature.append(result);arrays.update({f'q{order}_{k}':v for k,v in data.items()})
        path=destination/f'case_{i:02}.npz';np.savez_compressed(path,**arrays)
        result=dict(spec=record['spec'],original_full_plasma_status=record['original_full_plasma_status'],full_plasma_rate_accepted=False,
            input_arrays_sha256=record['arrays_sha256'],quadrature=quadrature,arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start)
        records.append(result);(destination/f'case_{i:02}.json').write_text(json.dumps(result,indent=2)+'\n')
        print('driving audit',record['spec'],result['original_full_plasma_status'],quadrature[-1]['metrics']['driving_driving_sum']['volume']['error_over_full_reference_rate_L2'],flush=True)
    if provenance(['examples.radial_drive_audit'])!=identity:raise ValueError('Source changed')
    (destination/'summary.json').write_text(json.dumps(dict(terminal=True,source_sha256=identity['source_sha256'],input_summary_sha256=sha(root/'summary.json'),cases=records),indent=2)+'\n')
