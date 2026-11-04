"""Bounded high-voltage startup gate, then actual microsecond HOPS and controls."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from examples.hybrid3d_campaign import freeze_source, write_json
from examples.iec_3d import DEFAULTS, run


def worker(args):
    out, root = Path(args.out), Path(args.root)
    state = dict(status='running', cases={}, scope='Collisionless 3D development; no optimized resonance or fusion result')
    write_json(out/'campaign.json', state)
    base = DEFAULTS | dict(bias_v=3000., drive_mode='dc', duration_ns=20., dt_ps=10.,
        particles=600, density_m3=1e11, source_ua=0., response_ns=2., saved_intervals=20)
    gates = []
    for dt in (10., 5.):
        name = 'startup_dt'+str(int(dt))
        config = base | dict(dt_ps=dt, run_label='HOPS high-voltage timestep comparison')
        target = out/name
        state['cases'][name] = dict(status='running', run=str(target))
        write_json(out/'campaign.json', state)
        result = run(target, config, args.mesh)
        h = json.loads((target/'history.json').read_text())
        particle_error = max(abs(r['energy_residual_j'])/max(r['electron_kinetic_j']+r['ion_kinetic_j']+r['absorbed_j']+r['injected_j'],1e-30) for r in h)
        gates.append(particle_error)
        state['cases'][name].update(status=result['status'], particle_energy_error=particle_error)
        write_json(out/'campaign.json', state)
        if result['status'] != 'completed':
            raise RuntimeError('High-voltage comparison did not complete')
    write_json(out/'startup_comparison.json', dict(dt_ps=[10,5], maximum_particle_energy_error=gates,
        note='Actual source-off +3 kV startups. Same initial sample; no particle pairing after absorption. Spatial force errors remain unresolved.'))
    if gates[1] > .05 or gates[1] > max(.005, 1.2*gates[0]):
        state.update(status='stopped_at_accuracy_gate', error='High-voltage timestep comparison needs investigation')
        write_json(out/'campaign.json', state)
        return
    long = base | dict(duration_ns=args.duration_us*1000, dt_ps=5., source_ua=1000.,
        response_ns=20., saved_intervals=500)
    for mode in ('hops', 'dc', 'fixed'):
        if (out/'STOP').exists():
            state['status'] = 'stopped'; break
        target = root/'outputs/lab_runs'/f'{out.name}_{mode}'
        config = long | dict(drive_mode=mode, run_label=f'HOPS · {mode} · +3 kV · 0–{args.duration_us:g} µs')
        state['cases'][mode] = dict(status='running', run=str(target), target_us=args.duration_us)
        write_json(out/'campaign.json', state)
        write_json(out/(mode+'_config.json'), config)
        result = run(target, config, args.mesh)
        state['cases'][mode].update(status=result['status'], reached_us=result['time_ns']/1000)
        write_json(out/'campaign.json', state)
        if result['status'] != 'completed':
            state['status'] = 'stopped'; break
    else:
        state['status'] = 'completed'
    write_json(out/'campaign.json', state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out', type=Path)
    parser.add_argument('--duration-us', type=float, default=5.)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--root', type=Path)
    args = parser.parse_args()
    if args.worker:
        try:
            return worker(args)
        except BaseException as error:
            p = args.out/'campaign.json'
            state = json.loads(p.read_text()) if p.exists() else {}
            for row in state.get('cases', {}).values():
                if row['status'] == 'running':
                    row.update(status='failed', error=str(error))
            state.update(status='failed', error=str(error))
            write_json(p, state)
            raise
    root = Path(__file__).resolve().parents[1]
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    mesh = args.mesh.resolve()
    shutil.copy2(mesh, args.out/'mesh.npz')
    mesh = args.out/'mesh.npz'
    source = freeze_source(root, args.out)
    command = [sys.executable, '-m', 'examples.hops_campaign', str(args.out), '--worker',
               '--root', str(root), '--mesh', str(mesh), '--duration-us', str(args.duration_us)]
    write_json(args.out/'launch.json', dict(command=command, source=str(source)))
    env = os.environ | dict(PYTHONPATH=str(source), OPENBLAS_NUM_THREADS='1',
                            OMP_NUM_THREADS='1', NUMBA_NUM_THREADS='1')
    with (args.out/'worker.log').open('w') as log:
        result = subprocess.run(command, cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError('HOPS worker failed; source and partial evidence retained')


if __name__ == '__main__':
    main()
