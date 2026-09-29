"""Replay raw fields, original equations and independently sampled local energies."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import numpy as np
from .variable_current_records import provenance,sha
from .variable_current_verification import initialize,pressure_data
from .variable_current_norms import measure
from plasma_bench.variable_compatible_current import weak_vector
from plasma_bench.kernel_induction import difference_curl


def relative(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),np.finfo(float).tiny))


def audit(directory):
    identity=json.loads((directory/'provenance.json').read_text())
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source identity differs')
    with tarfile.open(directory/'sources.tar.gz') as archive:
        for name,digest in identity['source_files'].items():
            if hashlib.sha256(archive.extractfile(name).read()).hexdigest()!=digest:raise RuntimeError('Archived source mismatch: '+name)
    records=json.loads((directory/'records.json').read_text());audits=[]
    for index,record in enumerate(records):
        if record['status']!='complete':raise RuntimeError('Unfinished case')
        path=directory/record['arrays_file']
        if sha(path)!=record['arrays_sha256']:raise RuntimeError('Raw array identity differs')
        with np.load(path) as data:arrays={k:data[k] for k in data.files}
        ref,m,field,nt,radial,axial=initialize(record['spec']);geo=m.geometry
        initial=dict(magnetic=field.packed_increment(geo),n_cell=m.density_map.cell_means,
                     nt_cell=nt,particle_radial=radial,particle_axial=axial,current=field.current(m))
        replay={key:relative(arrays[key]-value,value) for key,value in initial.items()}
        other=weak_vector(geo,ref.fields(geo.r,geo.z)['other_forcing'])
        density=m.density_forcing(arrays['current'],nt)
        replay['forcing_other']=relative(arrays['forcing_other']-other,other)
        replay['forcing_density']=relative(arrays['forcing_density']-density,density)
        pressure,diagnostics=pressure_data(m,field,ref,arrays['current'])
        replay.update({key:relative(arrays[key]-value,value) for key,value in pressure.items()})
        if max(replay.values())>1e-12:raise RuntimeError('Input or forcing replay mismatch')
        e=arrays['electric'];j=arrays['current'];jt=arrays['current_rate'];bt=arrays['magnetic_rate']
        inertia=m.inertia@jt;me=m.mass@e;f=other+density
        original=dict(ohm=relative(me-inertia-f,abs(me)+abs(inertia)+abs(f)),
            ampere=relative(m.mass@jt+geo.faraday.T@(m.hessian@bt),abs(m.mass@jt)+abs(geo.faraday.T@(m.hessian@bt))),
            faraday_backward=relative(difference_curl(geo,e)-bt,abs(geo.faraday)@abs(e)+abs(bt)),
            faraday_rate_relative=relative(difference_curl(geo,e)-bt,bt))
        if max(original[k] for k in ('ohm','ampere','faraday_backward'))>1e-10:raise RuntimeError('Original equation fails')
        norms,cells=measure(m.mesh,ref,geo.degree,arrays,order=14)
        independent={key:relative(arrays['production_'+key]-value,value) for key,value in cells.items()}
        if max(independent.values())>1e-7:raise RuntimeError('Independent physical local energy disagrees')
        density_pullback=float(np.sum(m.density_energy_gradient(j)*nt))
        density_rate=float(np.sum(m.mesh.volumes*cells['density_energy_rate']))
        pullback_scale=float(np.sum(m.mesh.volumes*abs(cells['density_energy_rate'])))
        pullback_error=abs(density_pullback-density_rate)/max(pullback_scale,np.finfo(float).tiny)
        if pullback_error>1e-7:raise RuntimeError('Density work pullback disagrees')
        result=dict(case=index,spec=record['spec'],input_replay=replay,original_equations=original,
            independent_energy_relative_difference=independent,density_pullback_backward=pullback_error,
            independent_norms=norms,pressure=diagnostics,arrays_sha256=sha(path))
        audits.append(result);print(index,original,flush=True)
    output=dict(source_sha256=identity['source_sha256'],cases=audits,
        records_sha256=sha(directory/'records.json'),source_archive_sha256=sha(directory/'sources.tar.gz'))
    (directory/'audit.json').write_text(json.dumps(output,indent=2)+'\n')
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    audit(parser.parse_args().directory)
