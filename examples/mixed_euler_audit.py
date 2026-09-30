"""Independent continuous equations, source work and mixed evolution replay."""
import argparse
import json
from pathlib import Path
import numpy as np
from .mixed_euler_records import provenance,sha
from .mixed_euler_study import initialize,totals
from .mixed_euler_step import step
from .mixed_euler_sampling import sample,grid,differences
from .cylindrical_entropy_reference import physical,rates


def audit(directory):
    identity=provenance();stored=json.loads((directory/'provenance.json').read_text())
    if stored['source_sha256']!=identity['source_sha256']:raise ValueError('Frozen source mismatch')
    results=[]
    for path in sorted(directory.glob('fluid_[0-9][0-9].json')):
        record=json.loads(path.read_text());raw=directory/record['arrays_file']
        if record['source_sha256']!=identity['source_sha256'] or sha(raw)!=record['arrays_sha256']:raise ValueError('Source or raw hash mismatch')
        a=dict(np.load(raw));ref,d,initial=initialize(record['spec']);u=a['state'];du,diag=d.derivative(u)
        if not np.array_equal(initial,a['initial']) or not np.array_equal(du,a['rate']):raise ValueError('Initial/RHS replay differs')
        n=record['accepted_steps'];dt=record['dt'];forced=record['spec']['manufactured'];h=a['history']
        if record['status']=='complete':
            if n!=record['steps']:raise ValueError('Incomplete completed run')
            replay,dw,ds=step(d,a['previous'],ref,(n-1)*dt,dt,forced)
            if not np.array_equal(replay,u) or a['work_before']+dw!=h[-1,8] or not np.array_equal(a['sources_before']+ds,h[-1,9:15]):raise ValueError('Last state/work/source replay differs')
        elif record['status']=='failed':
            before=u.copy()
            try:step(d,u,ref,n*dt,dt,forced)
            except ValueError as error:
                if repr(error)!=record['error']:raise ValueError('Different rejection')
            else:raise ValueError('Saved rejection did not recur')
            if not np.array_equal(before,u):raise ValueError('Rejected input mutated')
        else:raise ValueError('Nonterminal record')
        r,z,vol=grid(ref.radius,2*ref.half_length,32,256,24)
        rr,zz,vv=grid(ref.radius,2*ref.half_length,32,512,40);energies=[];angular=[];quadrature=[];alias=[]
        for state,index,row in zip(a['snapshots'],a['snapshot_steps'],h,strict=True):
            if not np.array_equal(totals(d,state,ref,index*dt,row[8],row[9:15]),row):raise ValueError('Physical history mismatch')
            energies_pair=[];angular_pair=[]
            for radius,points,volume in ((r,z,vol),(rr,zz,vv)):
                sampled=sample(d,state,radius,len(points));f=physical(sampled,radius[:,None],ref)
                if not all(np.isfinite(x).all() for x in f.values()) or np.min(f['pressure'])<=0:raise ValueError('Invalid independent physical sample')
                energies_pair.append(float(np.sum(volume*f['energy'])));angular_pair.append(float(np.sum(volume*sampled[...,2])))
            energies.append(energies_pair[-1]);angular.append(angular_pair[-1]);quadrature.append(abs(energies_pair[0]-energies_pair[1])/abs(energies_pair[1]));alias.append((row[1]-energies_pair[-1])/energies_pair[-1])
        low=sample(d,u,r,len(z));actual=physical(low,r[:,None],ref);total=float(np.sum(vol*actual['energy']))
        if total!=record['continuous_energy']:raise ValueError('Final physical energy differs')
        if forced:
            metrics=differences(low,ref.state(r[:,None],z[None,:],record['duration']),r,vol,ref,ref.radius/8)
            if metrics!=record['manufactured_errors']:raise ValueError('Manufactured physical norms differ')
        independent=rates(u,du,d.r,ref)['energy'];local=float(np.linalg.norm(independent-diag['physical_energy_rate'])/np.linalg.norm(independent))
        results.append(dict(file=path.name,status=record['status'],initial_rate_last_step_bitwise=True,
            independent_local_energy_rate_difference=local,maximum_quadrature_difference=max(quadrature),
            maximum_grid_continuous_energy_difference=max(abs(np.array(alias))),continuous_energy_history=energies,continuous_angular_history=angular,
            continuous_energy_change=(energies[-1]-energies[0])/energies[0],snapshots=len(h)))
        print(path.name,'audited',flush=True)
    result=dict(source_sha256=identity['source_sha256'],trajectories=results)
    (directory/'audit.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);args=parser.parse_args();audit(args.directory)


if __name__=='__main__':main()
