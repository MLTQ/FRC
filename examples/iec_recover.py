"""Reconstruct legacy restart metadata, requiring exact original source draws."""
import json
from pathlib import Path

import numpy as np

from plasma_bench.iec_checkpoint import digest, physics_hashes, write_checkpoint
from plasma_bench.iec_field import Field
from plasma_bench.iec_particles import Particles, WallEmitter, load_plasma
from plasma_bench.hops_drive import CageDrive


def exact(a,b):
    return a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes()


def recover(run, out):
    run,out=Path(run),Path(out)
    state=json.loads((run/'run.json').read_text()); config=state['config']
    history=json.loads((run/'history.json').read_text())
    provenance=json.loads((run/'provenance.json').read_text())
    hashes=physics_hashes()
    if provenance['numpy'] != np.__version__ or any(
        provenance['hashes']['plasma_bench/'+name+'.py']!=value for name,value in hashes.items()):
        raise ValueError('Legacy numerical library/physics differs; cannot claim exact reconstruction')
    field=Field(dict(np.load(run/'mesh.npz')))
    emitter=WallEmitter(field)
    rng=np.random.default_rng(config['seed'])
    p0=load_plasma(field,rng,config['particles'],config['seed_radius_m'],config['density_m3'],
                   config['electron_ev'],config['ion_ev'],config['ion_amu'])
    with np.load(run/'state_0000.npz') as f:
        if not all(exact(getattr(p0,k),f[k]) for k in ('x','v','cells','species')):
            raise RuntimeError('Reconstructed initial particle draw does not match')
    initial=float(p0.kinetic().sum()+history[0]['field_j'])
    steps=int(np.ceil(config['duration_ns']*1000/config['dt_ps']))
    dt=config['duration_ns']*1e-9/steps
    by_step={int(round(r['time_ns']*1e-9/dt)):i for i,r in enumerate(history)}
    last=max(by_step); previous=sorted(by_step)[-2]
    empty=Particles(np.empty((0,3)),np.empty((0,3)),np.empty(0,np.int64),
                    np.empty(0,np.int8),p0.weight,p0.ion_mass)
    emitted_charge,injected_energy,births=0.,0.,0
    checkpoints={}
    switches=json.loads((run/'drive_switches.json').read_text())
    for step in range(1,last+1):
        energy,charge=emitter.inject(empty,rng,config['source_ua']*1e-6,dt,config['source_ev'])
        injected_energy+=energy;emitted_charge+=charge;births+=len(empty.x)
        if len(empty.x): empty.retain(np.zeros(len(empty.x),bool))
        if step not in by_step: continue
        i=by_step[step]; row=history[i]
        expected=int(np.asarray(row['collected_macroparticles'])[:,0].sum())+row['electrons']-config['particles']
        if births!=expected or emitted_charge!=row['emitted_charge_c'] or injected_energy!=row['injected_j']:
            raise RuntimeError(f'Source-stream mismatch at original frame {i}')
        if step % 100000 == 0: print('source_reconstruction',step,step*dt*1e6,flush=True)
        if step not in (previous,last): continue
        with np.load(run/f'state_{i:04d}.npz') as f:
            p=Particles(*(f[k].copy() for k in ('x','v','cells','species')),float(f['macro_weight']),float(f['ion_mass']))
            phi,charge=f['phi'].copy(),f['charge'].copy()
            time_s=float(f['time_s'])
        q=field.deposit(p.x,p.cells,p.charge)
        solved,electric,electrodes,energy,space=field.solve(q,row['voltage_v'])
        if not exact(q,charge) or not exact(solved,phi) or energy!=row['field_j'] or space!=row['space_charge_j']:
            raise RuntimeError('Reconstructed Poisson state differs from original')
        drive=CageDrive(config);d=row['drive']
        drive.value=d['voltage_v'];drive.gate=d['gate_on']
        drive.filtered=d['filtered_radius_m'];drive.reference=d['reference_radius_m']
        drive.switches=[s for s in switches if s['time_us']<=time_s*1e6]
        event_time=lambda s:round(s['time_us']*1e-6/dt)*dt
        onsets=[event_time(s) for s in drive.switches if s['on']]
        drive.last_switch=event_time(drive.switches[-1]) if drive.switches else 0.
        drive.last_on=onsets[-1] if onsets else None
        drive.last_period=onsets[-1]-onsets[-2] if len(onsets)>1 else None
        drive.pulses=1+len(onsets)
        if drive.diagnostics()!=d: raise RuntimeError('Controller reconstruction mismatch')
        runtime=dict(p=p,rng=rng,drive=drive,history=history[:i+1],mesh_sha256=digest(run/'mesh.npz'),
            physics_hashes=hashes,phi=phi,charge=charge,electric=electric,electrodes=electrodes,
            lost_charge=np.asarray(row['collected_charge_c']),lost_count=np.asarray(row['collected_macroparticles'],dtype=np.int64),
            initial=initial,loss_energy=row['absorbed_j'],injected_energy=injected_energy,
            source_work=row['source_work_j'],injected_charge=emitted_charge,time_s=time_s,
            energy=energy,space=space,cage_voltage=row['voltage_v'])
        path=out/f'recovered_{step}.npz'
        write_checkpoint(path,config,step,dt,runtime);checkpoints[step]=path
    report=dict(initial_draw_exact=True,source_stream_frames_checked=len(history)-1,
                emitted_macroparticles=births,source_charge_exact=True,source_energy_exact=True,
                original_last_step=last,original_previous_step=previous,
                checkpoints={str(k):str(v) for k,v in checkpoints.items()})
    (out/'reconstruction.json').write_text(json.dumps(report,indent=2)+'\n')
    return config,previous,last,checkpoints
