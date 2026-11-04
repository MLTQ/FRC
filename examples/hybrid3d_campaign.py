"""Run a bounded, source-frozen queue of actual full-shot hybrid simulations."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time


CASES = (
    ('uniform_h2', 'Option A'), ('end_plugs', 'End plugs'),
    ('reboost', 'Re-boost'), ('inner_tube', 'Inner tube'),
    ('flux_core', 'Flux core'), ('core_punch', 'Core punch'),
)


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=1) + '\n')
    temporary.replace(path)


def freeze_source(project, out):
    """Freeze imports and rendering assets; output data remain outside git."""
    source = out / 'source'
    source.mkdir()
    for package in ('examples', 'plasma_bench'):
        shutil.copytree(project / package, source / package,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copytree(project / 'visualizations' / 'material_fields',
                    source / 'visualizations' / 'material_fields')
    hashes = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in source.rglob('*') if p.is_file()}
    write_json(out / 'source_hashes.json', hashes)
    diff = subprocess.run(['git', 'diff', '--binary'], cwd=project,
                          capture_output=True, check=True).stdout
    (out / 'source_diff.patch').write_bytes(diff)
    return source


def progress(run, row):
    """Read atomic diagnostics, including numerical sinks, without changing state."""
    history = run / 'history.json'
    if not history.exists():
        return None
    h = json.loads(history.read_text())
    first, last = h[0], h[-1]
    work = last.get('load_energy_work_j', 0.)
    error = last.get('load_energy_residual_j', 0.)
    row.update(reached_us=last['time_s'] * 1e6, steps=last['steps'],
               load_work_j=work, load_energy_residual_j=error,
               numerical_removed_energy_j=(last.get('cutoff_removed_heat_j', 0.)
                                            + last.get('electron_empty_node_energy_j', 0.)),
               loaded_ions_remaining_fraction=last['particles'] / first['particles'])
    if 'wall_s' in last:
        row['simulation_wall_s'] = last['wall_s']
    if abs(error) > max(5., .01 * abs(work)):
        return f'Load accounting error {error:g} J exceeds max(5 J, 1% of source work)'
    return None


def run_campaign(out, dx=.01, dz=.01, dt_ns=.5, workers=2, threads=2, gpu=0, gate_us=1.3, cases=None, fused_arrays=False):
    project = Path(__file__).resolve().parents[1]
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = freeze_source(project, out)
    prefix = out.name.replace('shots', 'full')
    env = dict(os.environ, NUMBA_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    rows = {}
    for case, label in CASES:
        if cases is not None and case not in cases:
            continue
        name = f'{prefix}_{case}'
        run = project / 'outputs' / name
        if run.exists():
            raise FileExistsError(run)
        backend_label = 'Metal' if gpu else 'CPU64'
        rows[name] = dict(status='queued', source=str(project / 'outputs' / f'fine_{case}'),
                          run=str(run), label=f'3D hybrid: {label} · 10 mm {backend_label} · 0–40 µs',
                          target_us=40., dx_m=dx, dz_m=dz, dt_ns=dt_ns,
                          ions_per_cell=8, floor_fraction=.05,
                          gpu_mode=gpu, electron_transport_cfl=.8, fused_cpu_arrays=fused_arrays,
                          scientific_status='development; grid and numerical-loss convergence unresolved')
    if not rows:
        raise ValueError('Select at least one known configuration')
    write_json(out / 'coordinator.json', dict(pid=os.getpid(), workers=workers,
               threads_per_worker=threads, gpu_mode=gpu, fused_cpu_arrays=fused_arrays, startup_gate_us=gate_us,
               source=str(source), started_unix=time.time(),
               initial_free_disk_gb=shutil.disk_usage(project).free / 1e9))
    active = {}
    try:
        while any(r['status'] in ('queued', 'running') for r in rows.values()):
            for name, (child, log) in list(active.items()):
                row = rows[name]
                problem = progress(Path(row['run']), row)
                if problem:
                    child.terminate()
                    child.wait(timeout=60)
                    row.update(status='invalidated', failure=problem)
                    write_json(Path(row['run']) / 'INVALIDATED.json', dict(reason=problem))
                code = child.poll()
                if code is None:
                    continue
                row['returncode'] = code
                log.close()
                del active[name]
                if row['status'] == 'invalidated':
                    continue
                if code != 0 or row.get('reached_us', 0.) < 40. - 1e-7:
                    row.update(status='failed', failure='Inspect simulation log; full shot incomplete')
                    continue
                row['status'] = 'complete_exploratory'
                # The exporter locks the shared run list and retains coarse runs.
                with (out / f'{name}_viewer.log').open('w') as viewer_log:
                    exported = subprocess.run([sys.executable, '-m', 'examples.hybrid3d_display',
                        str(project / 'outputs' / 'formation_display'),
                        row['run'] + '::' + row['label']], cwd=source, env=env,
                        stdout=viewer_log, stderr=subprocess.STDOUT)
                row['viewer_export_returncode'] = exported.returncode
                if exported.returncode:
                    row['status'] = 'export_failed'
                print(name, row['status'], flush=True)
            # Observe the corrected baseline past the previous1.141us failure
            # before filling the queue; do not discard this physical progress.
            first = next(iter(rows.values()))
            gate = first.get('reached_us', 0.) >= gate_us-1e-8
            if first['status'] in ('failed','invalidated'):
                for row in rows.values():
                    if row['status']=='queued':
                        row.update(status='deferred_baseline_failure', failure='Inspect baseline before starting more shots')
            limit = workers if gate else 1
            free = shutil.disk_usage(project).free
            if free < 10_000_000_000:
                for row in rows.values():
                    if row['status'] == 'queued':
                        row.update(status='deferred_disk', failure='Less than10GB free; no data deleted')
            for name, row in rows.items():
                if len(active) >= limit:
                    break
                if row['status'] != 'queued':
                    continue
                command = [sys.executable, '-m', 'examples.hybrid3d_frc', row['run'], row['source'],
                    't_us=0', f'dx={dx}', f'dz={dz}', 'ppc=8', 't_end_us=40',
                    'floor=0.05', 'drive=2', 'thermal=1', f'dt_ns={dt_ns}',
                    'snap_us=0.1', 'every_us=0.1', f'gpu={gpu}', f'arrays={int(fused_arrays)}']
                log = (out / f'{name}.log').open('w')
                child = subprocess.Popen(command, cwd=source, env=env,
                    stdout=log, stderr=subprocess.STDOUT)
                active[name] = (child, log)
                row.update(status='running', pid=child.pid, command=command,
                           started_unix=time.time())
                print('Started', name, child.pid, flush=True)
            write_json(out / 'campaign.json', rows)
            time.sleep(10)
        (out / 'DONE').write_text('Queue terminal; inspect campaign.json for failures.\n')
        return rows
    finally:
        # Only children created here are affected; never signal other simulations.
        for name, (child, log) in active.items():
            if child.poll() is None:
                child.terminate()
                try:
                    child.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
            log.close()
            rows[name]['status']='stopped_incomplete'
        if active:
            for row in rows.values():
                if row['status']=='queued':row['status']='deferred_interrupted'
            write_json(out/'campaign.json',rows)
            (out/'DONE').write_text('Stopped incomplete; inspect campaign.json.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out')
    parser.add_argument('--workers', type=int, choices=(1, 2), default=2)
    parser.add_argument('--threads', type=int, choices=(1, 2, 3, 4), default=2)
    parser.add_argument('--gpu', type=int, choices=(0,1,2), default=0)
    parser.add_argument('--arrays', action='store_true', help='Use the verified float64 fused energy stencils')
    parser.add_argument('--gate-us', type=float, default=1.3)
    parser.add_argument('--cases',nargs='+',choices=[name for name,_ in CASES])
    args = parser.parse_args()
    def terminate_owned_queue(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM,terminate_owned_queue)
    run_campaign(args.out, workers=args.workers, threads=args.threads, gpu=args.gpu,
                 gate_us=args.gate_us,cases=args.cases,fused_arrays=args.arrays)
