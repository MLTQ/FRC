"""Physical replays and whole-step timings for float64 electron loop fusion."""
import argparse
import cProfile
import hashlib
import json
import os
from pathlib import Path
import time

if __name__ == '__main__':
    for name in ('NUMBA_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
        os.environ.setdefault(name, '2' if name == 'NUMBA_NUM_THREADS' else '1')
import numpy as np
from examples.hybrid3d_frc import build, full_flux, diagnostics, save_snapshot
from examples.hybrid3d_heating_check import decay, entropy_contact, energy_advection
import plasma_bench.hybrid3d_energy_electrons as closure_module


def save(out, name, data):
    (out/f'{name}.json').write_text(json.dumps(data, indent=1)+'\n')


def independent(out, fused=False, gpu=0):
    options = dict(fused_arrays=True) if fused else {}
    backend = 'metal' if gpu else 'cpu'
    rows = []
    for cells in (32, 64, 128):
        rows.append(dict(kind='advection', **energy_advection(cells, backend=backend, **options)))
        save(out, 'independent', rows); print(rows[-1], flush=True)
    for cells in (16, 32, 64):
        rows.append(dict(kind='contact', **entropy_contact(cells, backend=backend,
                         gpu_moments=gpu >= 2, **options)))
        for eta in ('constant', 'spitzer'):
            rows.append(dict(kind=eta, **decay(cells, eta, backend=backend,
                             gpu_moments=gpu >= 2, **options)))
        save(out, 'independent', rows); print(rows[-3:], flush=True)


def device(out, fused=False, gpu=0, dx=.02, dz=.04, steps=1500, profile=False):
    project = Path(__file__).resolve().parents[1]
    state = full_flux(project/'outputs'/'fine_uniform_h2', 0.)
    options = dict(fused_arrays=True) if fused else {}
    start = time.perf_counter()
    model, info = build(state, dx, dz, ppc=8, floor=.05, seed=0, drive=2,
                        duration_s=40e-6, dt_ns=.5, thermal=True, gpu=gpu, **options)
    initialization = time.perf_counter()-start
    save(out, 'setup', info)
    history, measured = [diagnostics(model)], []
    profiler = cProfile.Profile() if profile else None
    for step in range(steps):
        if step == 10 and profiler is not None:
            profiler.enable()
        start = time.perf_counter()
        model.step(); model.source.advance(model)
        if step >= 10:
            measured.append(time.perf_counter()-start)
        if step == 9:
            history.append(diagnostics(model))
    if profiler is not None:
        profiler.disable(); profiler.dump_stats(str(out/'steps.prof'))
    history.append(diagnostics(model)); save(out, 'history', history)
    save_snapshot(model, out, 0, np.random.default_rng(1))
    np.savez_compressed(out/'state64.npz', b=np.asarray(model.b), x=model.x, v=model.v,
                        w=model.w, n=model.n, u=model.electrons.u,
                        flux_prev=model.flux_prev, flux_half=model.flux_half)
    result = dict(fused_arrays=fused, gpu=gpu, steps=steps, measured_steps=len(measured),
                  measured_s=sum(measured), initialization_s=initialization,
                  time_us=model.time*1e6, endpoint=history[-1])
    save(out, 'device', result); print(result, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out')
    parser.add_argument('--section', choices=('all','independent','device'), default='all')
    parser.add_argument('--fused', action='store_true')
    parser.add_argument('--gpu', type=int, choices=(0,1,2), default=0)
    parser.add_argument('--dx', type=float, default=.02)
    parser.add_argument('--dz', type=float, default=.04)
    parser.add_argument('--steps', type=int, default=1500)
    parser.add_argument('--profile', action='store_true')
    args = parser.parse_args()
    out = Path(args.out); out.mkdir(exist_ok=False, parents=True)
    package = Path(closure_module.__file__).parent
    files = [p for p in (package/'hybrid3d_arrays.py', package/'hybrid3d_energy_electrons.py',
             package/'hybrid3d_device.py', Path(__file__)) if p.exists()]
    save(out, 'source_hashes', {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    save(out, 'mode', dict(fused_arrays=args.fused, gpu=args.gpu,
                          closure_path=str(closure_module.__file__)))
    if args.section in ('all','independent'):
        independent(out, args.fused, args.gpu)
    if args.section in ('all','device'):
        device(out, args.fused, args.gpu, args.dx, args.dz, args.steps, args.profile)
