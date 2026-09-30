"""Replay actual rejected stages and audit their retained last accepted states."""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import numpy as np
from plasma_bench.log_fluid import LogFluidState
from .log_fluid_norms import FluidSampler


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    parser.add_argument('--radial-work',action='store_true');args=parser.parse_args();d=args.directory
    prefix='radial_work' if args.radial_work else 'log_fluid'
    study=importlib.import_module('examples.'+prefix+'_study');records=importlib.import_module('examples.'+prefix+'_records')
    identity=records.provenance();audits=[]
    for path in sorted(d.glob('fluid_*.json')):
        r=json.loads(path.read_text())
        if r['status']!='failed':continue
        assert r['source_sha256']==identity['source_sha256'];array_path=d/r['arrays_file']
        assert records.sha(array_path)==r['arrays_sha256'];a=np.load(array_path);ref,solver,_=study.initialize(r['spec'])
        state=LogFluidState(*(a[k+'_final'].copy() for k in study.STATE_KEYS))
        before=[getattr(state,k).copy() for k in study.STATE_KEYS]
        model,_,thermal,_,energy=solver.evaluate(state)
        history=np.array(study.totals(solver,state,r['last_accepted_time'],a['spatial_ledger_final']))
        mismatch=float(np.max(abs(history[:-2]-a['history'][-1,:-2])/np.maximum(abs(a['history'][-1,:-2]),1e-30)))
        assert mismatch<1e-13
        caught=None
        try:solver.step(state,r['last_accepted_time'],r['dt'])
        except (ValueError,RuntimeError) as error:caught=repr(error)
        assert caught==r['error'],(caught,r['error'])
        assert all(np.array_equal(v,getattr(state,k)) for v,k in zip(before,study.STATE_KEYS))
        quadrature=[];energies=[]
        for order in (18,22):
            grid=FluidSampler(solver.mesh,r['spec'].get('degree',5),order)
            minimum=float(grid.density(state.density).min());item=dict(order=order,min_density=minimum)
            if minimum<=0:
                item.update(positive_density=False,physical_energy_relative=None);energies.append(None)
            else:
                norms,independent=grid.measure(state,ref,r['last_accepted_time'])
                finite=bool(np.isfinite(independent).all());energies.append(independent if finite else None)
                error=float(np.linalg.norm(independent-energy)/np.linalg.norm(independent)) if finite else None
                item.update(positive_density=True,finite_energy=finite,physical_energy_relative=error,norms=norms)
            quadrature.append(item)
        convergence=(float(np.linalg.norm(energies[0]-energies[1])/np.linalg.norm(energies[1]))
                     if all(v is not None for v in energies) else None)
        warning=any(v['physical_energy_relative'] is None or v['physical_energy_relative']>2e-11 for v in quadrature)
        audits.append(dict(file=path.name,rejected_step_replayed=True,accepted_state_unchanged=True,
            history_relative=mismatch,last_state_quadrature=quadrature,higher_quadrature_relative_change=convergence,
            quadrature_warning=warning,
            min_density=float(model.density.min()),min_pressures=[float(f['pressure'].min()) for f in thermal]))
        print(path.name,'rejection/history verified; last-state quadrature warning',warning,flush=True)
    source=Path(__file__)
    (d/'failure_audit.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],audit_sha256=records.sha(source),records=audits),indent=2)+'\n')
    (d/'failure_audit_source.py').write_bytes(source.read_bytes());(d/'failure_audit_source.md').write_bytes(source.with_suffix('.md').read_bytes())
    print(len(audits),'actual failed steps replayed and last accepted states independently checked')


if __name__=='__main__':main()
