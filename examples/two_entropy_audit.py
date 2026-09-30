"""Exact replay and independent continuous physical integration of joint transport."""
import argparse
import json
from pathlib import Path
import numpy as np
from .two_entropy_records import provenance,sha
from .two_entropy_study import calculate_rate,initialize,totals
from .two_entropy_reference import measure


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    identity=provenance()['source_sha256'];assert json.loads((d/'provenance.json').read_text())['source_sha256']==identity
    rates=[];trajectories=[]
    for r in json.loads((d/'rates.json').read_text()):
        path=d/r['arrays_file'];assert sha(path)==r['arrays_sha256'] and r['source_sha256']==identity
        with np.load(path) as a:
            _,replay=calculate_rate(r['spec']);assert set(a.files)==set(replay) and all(np.array_equal(a[k],v) for k,v in replay.items())
        rates.append(dict(file=path.name,all_arrays_bitwise=True))
    for r in json.loads((d/'trajectories.json').read_text()):
        path=d/r['arrays_file'];assert sha(path)==r['arrays_sha256'] and r['source_sha256']==identity
        ref,solver,initial=initialize(r['spec'])
        with np.load(path) as a:
            assert np.array_equal(initial,a['initial']);assert np.array_equal(solver.derivative(initial)[0],a['initial_rate'])
            if r['accepted_steps']:assert np.array_equal(solver.step(a['before_final'],r['dt']),a['final'])
            if r['status']=='failed':
                before=a['final'].copy();caught=None
                try:solver.step(before,r['dt'])
                except (RuntimeError,ValueError) as error:caught=repr(error)
                assert caught==r['error'] and np.array_equal(before,a['final'])
            norms=[];quadrature=[]
            for i,(u,t) in enumerate(zip(a['snapshots'],a['times'])):
                h=totals(u,ref,solver,t);saved=a['history'][int(round(t/r['dt']))];assert np.array_equal(h,saved)
                high=measure(u,ref,t,16);low=measure(u,ref,t,8);norms.append(high)
                quadrature.append(abs(high['continuous_energy']-low['continuous_energy'])/abs(high['continuous_energy']))
            assert max(quadrature)<2e-12
        trajectories.append(dict(file=path.name,status=r['status'],initial_rate_last_step_history_bitwise=True,
                                 max_8x_vs_16x_energy_relative=max(quadrature),independent_16x_norms=norms))
        print(path.name,'audited',r['status'],flush=True)
    (d/'audit.json').write_text(json.dumps(dict(source_sha256=identity,rates=rates,trajectories=trajectories),indent=2)+'\n')


if __name__=='__main__':main()
