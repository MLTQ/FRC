"""Raw-source replay, original equations and independent physical thermal audits."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import numpy as np
from plasma_bench.kernel_induction import difference_curl
from .log_thermal_records import provenance,sha
from .log_thermal_cases import initialize
from .log_thermal_norms import measure
from .log_entropy_evolution import initialize as transport_initialize,physical_norms


def relative(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),np.finfo(float).tiny))


def verify_source(directory):
    identity=json.loads((directory/'provenance.json').read_text())
    if identity['source_sha256']!=provenance()['source_sha256']:raise RuntimeError('Live source differs from frozen study')
    with tarfile.open(directory/'sources.tar.gz') as archive:
        for name,digest in identity['source_files'].items():
            if hashlib.sha256(archive.extractfile(name).read()).hexdigest()!=digest:raise RuntimeError('Archived source differs: '+name)
    return identity


def load_arrays(directory,record):
    if record['status']!='complete':raise RuntimeError('Cannot audit an unfinished case')
    path=directory/record['arrays_file']
    if sha(path)!=record['arrays_sha256']:raise RuntimeError('Raw array identity differs')
    with np.load(path) as values:return {key:values[key] for key in values.files}


def audit_rates(directory):
    identity=verify_source(directory);records=json.loads((directory/'rates.json').read_text());results=[]
    for index,record in enumerate(records):
        a=load_arrays(directory,record);ref,m,field,nt,radial,axial,species,states=initialize(record['spec']);geo=m.geometry
        inputs=dict(n_cell=m.density_map.cell_means,nt_cell=nt,magnetic=field.packed_increment(geo),
            current=field.current(m),particle_radial=radial,particle_axial=axial)
        replay={key:relative(a[key]-value,value) for key,value in inputs.items()}
        for label,t in species.items():
            pair=t.pair(a['current'],states[label])
            for key,value in (('entropy',states[label]),('rate',pair.entropy_rate),('electric',pair.electric_full),
                              ('reduced',pair.electric_reduced),('potential',pair.potential),('energy_per_particle',pair.energy_per_particle)):
                replay[label+'_'+key]=relative(a[label+'_'+key]-value,value)
        if max(replay.values())>1e-12:raise RuntimeError('State or thermal rate replay differs')
        e=a['joint_electric'];jt=a['joint_current_rate'];bt=a['joint_magnetic_rate'];f=a['joint_forcing_full']
        me=m.mass@e;ij=m.inertia@jt;weak=geo.faraday.T@(m.hessian@bt)
        original=dict(ohm=relative(me-ij-f,abs(me)+abs(ij)+abs(f)),
            ampere=relative(m.mass@jt+weak,abs(m.mass@jt)+abs(weak)),
            faraday_backward=relative(difference_curl(geo,e)-bt,abs(geo.faraday)@abs(e)+abs(bt)),
            faraday_rate_relative=relative(difference_curl(geo,e)-bt,bt) if np.any(bt) else None,
            faraday_absolute=float(np.linalg.norm(difference_curl(geo,e)-bt)))
        if max(original[k] for k in ('ohm','ampere','faraday_backward'))>1e-10:raise RuntimeError('Original full induction equation fails')
        norms,cells=measure(m.mesh,ref,record['spec'],a,order=14)
        energy={label:{key:relative(a['production_'+label+'_'+key]-value,value) for key,value in values.items()} for label,values in cells.items()}
        if max(v for values in energy.values() for v in values.values())>1e-7:raise RuntimeError('Independent local energy/rate differs')
        results.append(dict(index=index,spec=record['spec'],input_replay=replay,original_equations=original,
            independent_energy=energy,independent_norms=norms,arrays_sha256=record['arrays_sha256']))
        print('rate audit',index,original,flush=True)
    output=dict(source_sha256=identity['source_sha256'],records_sha256=sha(directory/'rates.json'),cases=results)
    (directory/'audit.json').write_text(json.dumps(output,indent=2)+'\n');return output


def audit_transports(directory):
    identity=verify_source(directory);records=json.loads((directory/'transports.json').read_text());results=[]
    for index,record in enumerate(records):
        a=load_arrays(directory,record);ref,m,t,transport,s,magnetic=transport_initialize(record['spec'])
        for key,value in (('initial',s),('n_cell',m.density_map.cell_means),('current',transport.current),('magnetic',magnetic)):
            if not np.array_equal(a[key],value):raise RuntimeError('Transport initialization differs')
        if len(record['history'])!=record['steps']+1 or not np.isclose(record['steps']*record['dt'],ref.transit_time,rtol=2e-15):
            raise RuntimeError('Incomplete transport time history')
        last=transport.step(a['before_final'],record['dt'])
        if not np.array_equal(last.entropy,a['final']):raise RuntimeError('Last accepted RK4 step does not replay')
        if not np.array_equal(a['snapshots'][0],a['initial']) or not np.array_equal(a['snapshots'][-1],a['final']):
            raise RuntimeError('Transport snapshots do not span accepted endpoints')
        norms=physical_norms(m.mesh,ref,record['spec'],a,order=14)
        discrepancies=[]
        for physical in norms:
            step=int(round(physical['time']/record['dt']));stored=record['history'][step]
            discrepancies.append({key:abs(value-stored[key])/max(abs(value),np.finfo(float).tiny) for key,value in physical['integrals'].items()})
        if max(v for row in discrepancies for v in row.values())>1e-7:raise RuntimeError('Physical trajectory integral audit differs')
        if any(h['min_pressure']<=0 or not np.isfinite(list(h.values())).all() for h in record['history']):raise RuntimeError('Invalid accepted thermal history')
        results.append(dict(index=index,spec=record['spec'],independent_norms=norms,
            independent_integral_agreement=discrepancies,last_step_bitwise_replayed=True,
            all_saved_snapshots_audited=True,arrays_sha256=record['arrays_sha256']))
        print('transport audit',index,flush=True)
    output=dict(source_sha256=identity['source_sha256'],records_sha256=sha(directory/'transports.json'),cases=results)
    (directory/'audit.json').write_text(json.dumps(output,indent=2)+'\n');return output


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('rates',type=Path);parser.add_argument('transports',type=Path)
    args=parser.parse_args();audit_rates(args.rates);audit_transports(args.transports)
