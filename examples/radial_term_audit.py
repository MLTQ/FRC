"""Independent component and radial-region errors of the saved full RHS."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.strong_wall_plasma import StrongWallFluid
from plasma_bench.mhd_mesh import MHDMesh
from .archived_contact_mode import ArchivedContactMode
from .radial_parent_terms import terms
from .current_transport_sampling import TransportSampler
from .radial_term_sampling import radial_slice,velocity_rate,negative_curl,raw_faraday,projected_faraday,wall_projection
from .implicit_plasma_records import provenance,freeze,sha


ELECTRIC=('ideal','hall','density','transport','pressure_full','inertia')
CURL=ELECTRIC+('pressure_reduced','pressure_potential')
FLUID=('bulk','lorentz','drift','species','wall_reaction')


def measure(record,raw,order):
    mode=ArchivedContactMode();ref=mode.ref;spec=record['spec'];epsilon=record['epsilon']
    d=StrongWallFluid(ref.radius,ref.length,spec['elements'],spec['points'],ref.ion_mass+ref.electron_mass,
        degree=5,gammas=(ref.gamma_e,ref.gamma_i));mesh=MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    sample=TransportSampler(d,mesh,order);projection=wall_projection(mesh,spec.get('quadrature',12))
    center=lambda key:(raw['plus_'+key]-raw['minus_'+key])/(2*epsilon)
    ecoeff={k:center('electric_'+k) for k in CURL};ecoeff['total']=center('electric')
    bcoeff={k:projected_faraday(mesh,v,projection) for k,v in ecoeff.items()}
    bcoeff['sum_full']=sum(bcoeff[k] for k in ELECTRIC)
    bcoeff['sum_reduced']=sum(bcoeff[k] for k in ELECTRIC if k!='pressure_full')+bcoeff['pressure_reduced']
    bcoeff['accepted_rate']=center('magnetic_rate')
    profiles={};hashes={};cross_profiles={};net_profiles={}
    for start in range(0,mesh.nr,8):
        stop=min(start+8,mesh.nr);s=radial_slice(sample,start*order,stop*order);count=stop-start
        reference=terms(mode,s.R,s.Z);target=mode.fields(s.R,s.Z)
        nets=dict(momentum=mode.growth*target['physical'][...,1:4],electric=target['electric'],faraday=mode.growth*target['physical'][...,6:9])
        def integral(v):return np.sum((s.weight[...,None]*v).reshape(count,order,len(s.z),3),axis=(1,2))
        errors={}
        def add(group,name,actual,expected):
            error=actual-expected;key=group+'_'+name
            values=np.stack((integral(actual**2),integral(expected**2),integral(error**2)),axis=1)
            if key not in profiles:profiles[key]=np.zeros((mesh.nr,3,3));hashes[key]=hashlib.sha256()
            profiles[key][start:stop]=values;hashes[key].update(np.ascontiguousarray(actual).tobytes())
            errors[key]=error
        for group,net in nets.items():
            if group not in net_profiles:net_profiles[group]=np.zeros((mesh.nr,3))
            net_profiles[group][start:stop]=integral(net**2)
        momentum={}
        for name in FLUID:
            value=(velocity_rate(s,raw['plus_state'],raw['plus_fluid_'+name])-velocity_rate(s,raw['minus_state'],raw['minus_fluid_'+name]))/(2*epsilon)
            momentum[name]=value;expected=reference['momentum'].get('pressure' if name=='bulk' else name,np.zeros_like(value))
            add('momentum',name,value,expected)
        add('momentum','sum',sum(momentum.values()),nets['momentum'])
        accepted=(velocity_rate(s,raw['plus_state'],raw['plus_fluid_rate'])-velocity_rate(s,raw['minus_state'],raw['minus_fluid_rate']))/(2*epsilon)
        add('momentum','accepted_rate',accepted,nets['momentum'])
        for name in ELECTRIC+('total',):
            value=s.electric(ecoeff[name]);expected=(nets['electric'] if name=='total' else reference['electric'].get('pressure' if name=='pressure_full' else name,np.zeros_like(value)))
            add('electric',name,value,expected)
        for name,coefficient in bcoeff.items():
            value=s.magnetic(coefficient)
            if name in ('total','sum_full','sum_reduced','accepted_rate'):expected=nets['faraday']
            else:expected=reference['induction'].get('pressure' if name.startswith('pressure_') and name!='pressure_potential' else name,np.zeros_like(value))
            add('faraday',name,value,expected)
        for name in ('ideal','hall','pressure_full','pressure_reduced','inertia'):
            value=negative_curl(s,ecoeff[name]);expected=reference['induction']['pressure' if name.startswith('pressure_') else name]
            add('strong_faraday',name,value,expected)
        for group,names in [('momentum',FLUID),('faraday',('ideal','hall','pressure_reduced','inertia','density','transport'))]:
            key=group+'_error_cross';count_terms=len(names)
            if key not in cross_profiles:cross_profiles[key]=np.zeros((mesh.nr,count_terms,count_terms,3))
            for i,a in enumerate(names):
                for j,b in enumerate(names):cross_profiles[key][start:stop,i,j]=integral(errors[group+'_'+a]*errors[group+'_'+b])
    metrics={}
    for key,values in profiles.items():
        group=key.split('_')[0];group='faraday' if group=='strong' else group;net=net_profiles[group]
        regions={}
        for name,indices in [('volume',slice(None)),('axis',slice(0,mesh.nr//d.elements)),('wall',slice(mesh.nr-mesh.nr//d.elements,None))]:
            a,t,e=np.sum(values[indices],axis=0);n=np.sum(net[indices],axis=0)
            regions[name]=dict(actual_L2=np.sqrt(a).tolist(),target_L2=np.sqrt(t).tolist(),error_L2=np.sqrt(e).tolist(),
                error_over_net_L2=np.sqrt(e/np.maximum(n,np.finfo(float).tiny)).tolist(),
                relative_to_own_term=[float(np.sqrt(x/y)) if y>0 else None for x,y in zip(e,t)])
        metrics[key]=regions
    arrays={**profiles,**cross_profiles,**{'net_'+k:v for k,v in net_profiles.items()},'radius':mesh.r}
    return dict(order=order,metrics=metrics,physical_stream_sha256={k:v.hexdigest() for k,v in hashes.items()},
        cross_term_order=dict(momentum=list(FLUID),faraday=['ideal','hall','pressure_reduced','inertia','density','transport'])),arrays


def run(destination,original='outputs/radial16_term_capture'):
    destination=Path(destination);root=Path(original)
    if destination.exists():raise ValueError('Refusing overwrite')
    source=json.loads((root/'summary.json').read_text())
    if not source['terminal']:raise ValueError('Terminal captured input required')
    identity=provenance(['examples.radial_term_audit']);destination.mkdir(parents=True);freeze(destination,identity);records=[]
    for index,record in enumerate(source['cases']):
        start=time.monotonic();path=root/record['arrays_file']
        if sha(path)!=record['arrays_sha256']:raise ValueError('Captured input changed')
        raw=dict(np.load(path));quadrature=[];arrays={}
        for order in (16,20):
            measured,data=measure(record,raw,order);quadrature.append(measured);arrays.update({f'q{order}_{k}':v for k,v in data.items()})
        path=destination/f'case_{index:02}.npz';np.savez_compressed(path,**arrays)
        result=dict(spec=record['spec'],input_arrays_sha256=record['arrays_sha256'],quadrature=quadrature,arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start)
        records.append(result);(destination/f'case_{index:02}.json').write_text(json.dumps(result,indent=2)+'\n')
        print('parent-term audit',record['spec'],quadrature[-1]['metrics']['faraday_accepted_rate']['volume']['error_over_net_L2'],flush=True)
    if provenance(['examples.radial_term_audit'])!=identity:raise ValueError('Source changed')
    (destination/'summary.json').write_text(json.dumps(dict(terminal=True,source_sha256=identity['source_sha256'],input_summary_sha256=sha(root/'summary.json'),cases=records),indent=2)+'\n')
