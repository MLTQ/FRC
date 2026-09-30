"""Independent physical and original-platform replay of regular-coordinate evidence."""
import argparse
import json
from pathlib import Path
import numpy as np
from .regular_euler_records import provenance,sha
from .regular_euler_study import initialize,totals,rate_case
from .regular_euler_spectrum import calculate
from .regular_euler_reference import energy
from .regular_euler_sampling import sample
from .regular_euler_step import step
from .mixed_euler_sampling import grid,differences
from .cylindrical_entropy_reference import physical


def audit(directory):
    identity=provenance()
    if json.loads((directory/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:raise ValueError('Frozen source mismatch')
    result=[]
    for path in sorted(directory.glob('*.json')):
        if path.stem in ('provenance','audit') or 'checkpoint' in path.stem:continue
        record=json.loads(path.read_text());raw=directory/record['arrays_file']
        if record['source_sha256']!=identity['source_sha256'] or sha(raw)!=record['arrays_sha256']:raise ValueError('Source/raw hash mismatch')
        a=dict(np.load(raw));entry=dict(file=path.name)
        if path.stem.startswith(('rates_','spectra_')):
            replay,arrays=rate_case(record['spec']) if path.stem.startswith('rates_') else calculate(record['elements'],record['degree'])
            if any(not np.array_equal(value,a[key]) for key,value in arrays.items()):raise ValueError('Non-bitwise rate/spectrum replay')
            if any(value!=record[key] for key,value in replay.items()):raise ValueError('Rate/spectrum metric mismatch')
            entry['replayed_arrays']=len(arrays)
        else:
            ref,d,initial=initialize(record['spec']);u=a['state'];du,diag=d.derivative(u);n=record['accepted_steps'];dt=record['dt'];h=a['history'];forced=record['spec']['manufactured']
            if not np.array_equal(initial,a['initial']) or not np.array_equal(du,a['rate']):raise ValueError('Initial/RHS mismatch')
            if record['status']=='complete':
                if n!=record['steps']:raise ValueError('Incomplete completed run')
                replay,dw,ds=step(d,a['previous'],ref,(n-1)*dt,dt,forced)
                if not np.array_equal(replay,u) or a['work_before']+dw!=h[-1,8] or not np.array_equal(a['sources_before']+ds,h[-1,9:15]):raise ValueError('Last state/work/source mismatch')
            elif record['status']=='failed':
                before=u.copy()
                try:step(d,u,ref,n*dt,dt,forced)
                except ValueError as error:
                    if repr(error)!=record['error']:raise ValueError('Different rejection')
                else:raise ValueError('Rejection did not recur')
                if not np.array_equal(u,before):raise ValueError('Rejected step mutated state')
            else:raise ValueError('Nonterminal record')
            r,z,v=grid(ref.radius,2*ref.half_length,32,256,24);rr,zz,vv=grid(ref.radius,2*ref.half_length,32,512,40)
            energies=[];angular=[];quadrature=[];alias=[];max_speed=[]
            for state,index,row in zip(a['snapshots'],a['snapshot_steps'],h,strict=True):
                if not np.array_equal(totals(d,state,ref,index*dt,row[8],row[9:15]),row):raise ValueError('Physical history mismatch')
                e=[];l=[]
                for radius,points,volume in ((r,z,v),(rr,zz,vv)):
                    field=sample(d,state,radius,len(points));f=physical(field,radius[:,None],ref)
                    if not all(np.isfinite(x).all() for x in f.values()) or np.min(f['pressure'])<=0:raise ValueError('Invalid independent physical state')
                    e.append(float(np.sum(volume*f['energy'])));l.append(float(np.sum(volume*field[...,2])))
                energies.append(e[-1]);angular.append(l[-1]);quadrature.append(abs(e[0]-e[1])/abs(e[1]));alias.append((row[1]-e[-1])/e[-1]);max_speed.append(float(np.max(abs(f['velocity']))/ref.sound))
            low=sample(d,u,r,len(z));f=physical(low,r[:,None],ref)
            if float(np.sum(v*f['energy']))!=record['continuous_energy']:raise ValueError('Final independent energy mismatch')
            if forced and differences(low,ref.state(r[:,None],z[None,:],record['duration']),r,v,ref,ref.radius/8)!=record['manufactured_errors']:raise ValueError('Manufactured errors mismatch')
            independent=np.imag(energy(u+1j*1e-28*du,d.x,ref))/1e-28
            entry.update(status=record['status'],initial_rate_last_step_bitwise=True,independent_local_energy_rate_difference=float(np.linalg.norm(independent-diag['physical_energy_rate'])/np.linalg.norm(independent)),
                continuous_energy_history=energies,continuous_angular_history=angular,maximum_quadrature_difference=max(quadrature),maximum_grid_continuous_energy_difference=max(abs(np.array(alias))),
                maximum_sampled_speed_over_sound=max(max_speed),continuous_energy_change=(energies[-1]-energies[0])/energies[0],snapshots=len(h))
        result.append(entry);print(path.name,'audited',flush=True)
    output=dict(source_sha256=identity['source_sha256'],records=result);(directory/'audit.json').write_text(json.dumps(output,indent=2)+'\n');return output


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);args=p.parse_args();audit(args.directory)


if __name__=='__main__':main()
