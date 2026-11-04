"""Bounded actual-device replays isolating Metal transport and CIC precision."""
import argparse
import gc
import json
import os
from pathlib import Path
import subprocess
import sys
import time

if __name__ == '__main__':
    os.environ.setdefault('NUMBA_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    os.environ.setdefault('OMP_NUM_THREADS', '1')

import numpy as np
from examples.hybrid3d_campaign import freeze_source, write_json
from examples.hybrid3d_frc import full_flux, build, diagnostics
from examples.hybrid3d_gpu_check import independent_cic
from plasma_bench.hybrid3d import deposit, faces_to_nodes
from plasma_bench.hybrid3d_metal import MetalKernels
from plasma_bench.hybrid3d_metal_cic import SplitCoordinateMetal


MODES = dict(cpu=(0, False), cpu_repeat=(0, False), transport=(1, False),
             deposition=(0, True), combined=(2, False),
             split_deposition=(0, 'split'), split_combined=(1, 'split'))


def difference(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    delta = a-b
    scale = float(np.linalg.norm(b.ravel()))
    return dict(relative_l2=float(np.linalg.norm(delta.ravel())/scale) if scale else None,
                absolute_max=float(np.max(np.abs(delta))) if delta.size else 0.,
                bit_identical=bool(np.array_equal(a, b)))


def moment_probe(model, out, name):
    """Actual native particles; expose casting and nonrepeatable accumulation."""
    gpu = MetalKernels(model)
    cpu = deposit(model.x, model.v, model.w, model.lo, model.d, model.shape, model.nthreads)
    quantized = deposit(*(np.ascontiguousarray(a, dtype=np.float32).astype(float)
                          for a in (model.x, model.v, model.w)),
                        model.lo, model.d, model.shape, model.nthreads)
    first = gpu.moments(model.x, model.v, model.w)
    repeats = [difference(gpu.moments(model.x, model.v, model.w), first) for _ in range(3)]
    expected_number = float(model.w.sum())
    expected_flux = np.sum(model.w[:, None]*model.v, axis=0)
    def measures(moments):
        return dict(density=difference(moments[0], cpu[0]), flux=difference(moments[1:], cpu[1:]),
            integrated_number_relative_error=float((moments[0].sum()*model.box.cell_volume-expected_number)/expected_number),
            integrated_flux_relative_error=float(np.linalg.norm(moments[1:].sum(axis=(1,2,3))*model.box.cell_volume-expected_flux)
                                               /max(np.linalg.norm(expected_flux), 1e-300)))
    write_json(out/f'moments_{name}.json', dict(ions=len(model.x), time_s=model.time,
               cpu=measures(cpu), gpu=measures(first), quantized_input_cpu=measures(quantized),
               gpu_repeat_differences=repeats,
               note='Input quantization retains CPU64 normalization/accumulation; it is a sensitivity decomposition, not the GPU arithmetic.'))


def snapshot(model, out, step):
    """Native float64 fields and full particles, not a viewer's float32 sample."""
    response = np.asarray(model.b)-model.source(model.time)
    np.savez_compressed(out/f'state_{step:05d}.npz', time_s=model.time,
        b=np.asarray(model.b), response=response, n=model.n, pe=model._pressure(model.n),
        u=model.electrons.u, x=model.x, v=model.v, w=model.w,
        flux_prev=model.flux_prev, flux_half=model.flux_half)


def compare(out, modes, reference=None):
    reference = Path(reference).resolve() if reference else out/'cpu'
    if not reference.exists():
        return
    info = json.loads((reference/'setup.json').read_text())
    archive = []
    if not list(reference.glob('state_*.npz')):
        for path in reference.glob('snap_*.npz'):
            with np.load(path) as frame:
                archive.append((float(frame['time_s']), path))
    rows = []
    for name in modes:
        if name == 'cpu' or not (out/name/'history.json').exists():
            continue
        for path in sorted((out/name).glob('state_*.npz')):
            other = reference/path.name
            archived = False
            if not other.exists():
                if not archive:
                    continue
                with np.load(path) as frame:
                    target = float(frame['time_s'])
                time_s, other = min(archive, key=lambda pair:abs(pair[0]-target))
                if abs(time_s-target) > 1e-12:
                    continue
                archived = True
            with np.load(path) as a, np.load(other) as raw:
                b = raw if not archived else dict(time_s=raw['time_s'],
                    b=np.stack([raw[key] for key in ('bx','by','bz')]).astype(float),
                    response=raw['b_response'].astype(float), n=raw['n'].astype(float),
                    pe=raw['pe'].astype(float),
                    u=raw['pe'].astype(float)/(info['thermal_electrons']['gamma']-1))
                if abs(float(a['time_s'])-float(b['time_s'])) > 1e-15:
                    raise RuntimeError('Unmatched comparison times')
                if a['n'].shape != b['n'].shape:
                    raise RuntimeError('Mismatched comparison grids')
                row = dict(mode=name, time_us=float(a['time_s'])*1e6,
                           reference_precision='saved_float32' if archived else 'native_float64',
                           fields={key:difference(a[key], b[key]) for key in
                                   ('b','response','n','pe','u','flux_prev','flux_half') if key in b})
                if archived:
                    row['particles'] = dict(note='Reference has only a display sample; no particle pairing')
                elif a['x'].shape == b['x'].shape:
                    row['particles'] = {key:difference(a[key], b[key]) for key in ('x','v','w')}
                else:
                    row['particles'] = dict(note='Absorption histories differ; no identity pairing asserted',
                                            count_a=len(a['x']),count_cpu=len(b['x']))
                mean = .5*(a['n']+b['n'])
                # The free-node support is copied once from the actual device.
                support = np.load(out/'node_gas_fraction.npy') > 0
                masks = dict(dense=support & (mean >= info['floor_m3']),
                             thin=support & (mean > 0) & (mean < info['floor_m3']),
                             empty=support & (mean == 0))
                aa, bb = faces_to_nodes(*a['response']), faces_to_nodes(*b['response'])
                row['node_regions'] = {key:dict(nodes=int(mask.sum()),
                    field=difference(aa[:, mask], bb[:, mask]), density=difference(a['n'][mask], b['n'][mask]))
                    for key,mask in masks.items() if mask.any()}
                rows.append(row)
    write_json(out/'comparison.json', dict(reference='cpu', rows=rows,
                note='CPU is a numerical reference, not physical truth; unweighted region norms are not energy integrals. Saved float32 references limit small-difference resolution and cannot establish bitwise agreement.'))


def worker(args):
    out = Path(args.out).resolve()
    state = full_flux(Path(args.input).resolve(), 0.)
    progress = dict(status='running', pid=os.getpid(), target_us=args.end_us, cases={})
    write_json(out/'progress.json', progress)
    try:
        for name in args.modes:
            run = out/name
            run.mkdir()
            started = time.monotonic()
            mode, moments_only = MODES[name]
            model, info = build(state, args.dx, args.dz, ppc=8, floor=.05, seed=0,
                drive=2, duration_s=40e-6, dt_ns=args.dt_ns, thermal=True, gpu=mode)
            if moments_only:
                cls = SplitCoordinateMetal if moments_only == 'split' else MetalKernels
                model._metal_moments = cls(model)
                info['compute_backend'] = ('split_coordinate_metal_cic' if moments_only == 'split' else 'metal_cic')
                info['thermal_electrons']['gpu_moments'] = True
                info['thermal_electrons']['backend'] = 'metal' if mode else 'cpu'
            info['precision_replay'] = dict(mode=name, dx=args.dx, dz=args.dz, seed=0)
            write_json(run/'setup.json', info)
            if not (out/'node_gas_fraction.npy').exists():
                np.save(out/'node_gas_fraction.npy', model.node_gas_fraction)
            steps = int(round(args.end_us*1e-6/model.dt))
            if abs(steps*model.dt-args.end_us*1e-6) > 1e-14:
                raise ValueError('Endpoint must be an integral number of ion steps')
            stride = max(1, steps//2)
            history = [diagnostics(model)]
            snapshot(model, run, 0)
            if name == 'cpu':
                np.save(out/'node_gas_fraction.npy', model.node_gas_fraction)
                moment_probe(model, out, 'initial')
                independent_cic(out)
            progress['cases'][name] = dict(status='running', time_us=0., steps=steps)
            write_json(out/'progress.json', progress)
            for step in range(1, steps+1):
                model.step()
                model.source.advance(model)
                if step % stride == 0 or step == steps:
                    row = diagnostics(model)
                    row['elapsed_s'] = time.monotonic()-started
                    history.append(row)
                    snapshot(model, run, step)
                    write_json(run/'history.json', history)
                    progress['cases'][name].update(time_us=model.time*1e6, reached_steps=step,
                        elapsed_s=row['elapsed_s'], source_work_j=row['load_energy_work_j'],
                        ledger_residual_j=row['load_energy_residual_j'])
                    write_json(out/'progress.json', progress)
                    print(name, step, model.time*1e6, row['elapsed_s'], flush=True)
                    if abs(row['load_energy_residual_j']) > max(5., .01*abs(row['load_energy_work_j'])):
                        raise RuntimeError('Load ledger exceeds the existing campaign accounting guard; state retained without correction')
            if name == 'cpu':
                moment_probe(model, out, 'endpoint')
            progress['cases'][name]['status'] = 'completed'
            write_json(out/'progress.json', progress)
            del model
            gc.collect()
            compare(out, args.modes, args.reference)
        progress['status'] = 'completed'
        write_json(out/'progress.json', progress)
    except BaseException as error:
        progress.update(status='failed', error=str(error))
        write_json(out/'progress.json', progress)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out')
    parser.add_argument('--input', required=True)
    parser.add_argument('--dx', type=float, default=.01)
    parser.add_argument('--dz', type=float, default=.01)
    parser.add_argument('--dt-ns', type=float, default=.5)
    parser.add_argument('--end-us', type=float, default=.1)
    parser.add_argument('--modes', nargs='+', choices=tuple(MODES),
                        default=['cpu','cpu_repeat','transport','deposition','combined'])
    parser.add_argument('--reference', help='Existing matched CPU run, for prototype-only additions')
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    project = Path(__file__).resolve().parents[1]
    source = freeze_source(project, out)
    command = [sys.executable, '-m', 'examples.hybrid3d_precision_check', str(out), '--worker',
               '--input', str(Path(args.input).resolve()), '--dx', str(args.dx), '--dz', str(args.dz),
               '--dt-ns', str(args.dt_ns), '--end-us', str(args.end_us), '--modes', *args.modes]
    if args.reference:
        command += ['--reference', str(Path(args.reference).resolve())]
    write_json(out/'launch.json', dict(command=command, source=str(source), pid=os.getpid(),
                                     scope='Bounded initial-fill numerical sensitivity; no full campaign'))
    environment = os.environ.copy()
    environment['PYTHONPATH'] = str(source)
    with (out/'worker.log').open('w') as log:
        result = subprocess.run(command, cwd=source, env=environment, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Bounded precision comparison failed; inspect {out}/worker.log')
    print(out/'comparison.json', flush=True)


if __name__ == '__main__':
    main()
