"""Archived baseline equivalence, amplitude and quadrature controls."""
import json
from pathlib import Path
import numpy as np
from .implicit_plasma_records import provenance,freeze,sha


def load(directory,index):
    root=Path('outputs')/directory;record=json.loads((root/'summary.json').read_text())['cases'][index]
    if record['status']!='complete':raise ValueError('Complete case required')
    path=root/record['arrays_file']
    if sha(path)!=record['arrays_sha256']:raise ValueError('Raw case hash mismatch')
    return record,dict(np.load(path)),str(path)


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    old,a,oldpath=load('radial16_matched_rates',5)
    base,b,basepath=load('radial16_balanced_rates',0)
    if a.keys()!=b.keys() or any(not np.array_equal(a[k],b[k]) for k in a):
        raise ValueError('Baseline changed')
    first,a,firstpath=load('radial16_balanced_rates',1);second,b,secondpath=load('radial16_balanced_rates',2)
    if any(first['spec'][k]!=second['spec'][k] for k in ('elements','points','field_nr','field_nz')):raise ValueError('Same grid required')
    arrays={};amp={}
    for key in ('fluid_rate','magnetic_rate','current_rate','electric'):
        x=(a['plus_'+key]-a['minus_'+key])/(2*first['epsilon'])
        y=(b['plus_'+key]-b['minus_'+key])/(2*second['epsilon'])
        arrays[key+'_difference']=x-y
        amp[key]=float(np.linalg.norm(x-y)/np.linalg.norm(y))
    quadrature=[]
    for name in ('radial16_matched_rates','radial16_fine_rates','radial16_balanced_rates'):
        data=json.loads((Path('outputs')/name/'summary.json').read_text())
        if not data['terminal']:raise ValueError('Terminal batch required')
        for index,case in enumerate(data['cases']):
            if case['status']!='complete':continue
            low,high=case['quadrature'];errors=[]
            for group in ('rate_errors','initial_errors'):
                for key in low[group]:
                    for quantity in ('relative_L2','axis_element','wall_element'):
                        x=low[group][key][quantity];y=high[group][key][quantity]
                        errors.append(abs(x-y)/max(abs(x),abs(y),np.finfo(float).tiny))
            quadrature.append(dict(directory=name,index=index,max_relative_physical_norm_difference=max(errors)))
    identity=provenance(['examples.radial16_comparison']);destination.mkdir(parents=True);freeze(destination,identity)
    path=destination/'amplitude.npz';np.savez_compressed(path,**arrays)
    result=dict(terminal=True,source_sha256=identity['source_sha256'],baseline_all_arrays_bitwise_equal=True,
        input_hashes={p:sha(p) for p in (oldpath,basepath,firstpath,secondpath)},amplitude_coefficient_relative_difference=amp,
        amplitude_metric='Euclidean coefficient norms; not physical reconstruction norms',quadrature=quadrature,
        arrays_file=path.name,arrays_sha256=sha(path))
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('comparison',amp,'quadrature max',max(r['max_relative_physical_norm_difference'] for r in quadrature),flush=True)
