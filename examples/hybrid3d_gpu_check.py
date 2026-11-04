"""Physical Metal-backend comparisons and synchronized whole-step timings."""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace
if __name__ == '__main__':
    os.environ.setdefault('NUMBA_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from examples.hybrid3d_heating_check import entropy_contact, decay, energy_advection
from examples.hybrid3d_frc import full_flux, build, diagnostics, save_snapshot
from plasma_bench.hybrid3d import Box, deposit
from plasma_bench.hybrid3d_metal import MetalKernels


def save(out, name, rows):
    (out/f'{name}.json').write_text(json.dumps(rows, indent=1)+'\n')


def independent(out):
    independent_cic(out)
    advection=[]
    for backend in ('cpu','metal'):
        for cells in (32,64,128):
            advection.append(energy_advection(cells,backend=backend))
            save(out,'energy_advection',advection);print('energy advection',advection[-1],flush=True)
    rows = []
    for gpu_moments in (False, True):
        for cells in (16, 32, 64):
            contact = entropy_contact(cells, backend='metal', gpu_moments=gpu_moments)
            rows.append(dict(kind='moving_contact', gpu_moments=gpu_moments, **contact))
            save(out, 'independent', rows); print(rows[-1], flush=True)
            for eta in ('constant', 'spitzer'):
                rows.append(dict(kind='force_free_heating', gpu_moments=gpu_moments,
                    **decay(cells, eta, backend='metal', gpu_moments=gpu_moments)))
                save(out, 'independent', rows); print(rows[-1], flush=True)


def independent_cic(out):
    """Uniform density/drift and exact integrated number, not a CPU-kernel oracle."""
    box = Box((16,16,16),(0.,0.,0.),(.25,.25,.25))
    density = 2e20
    grid = np.indices(box.shape).reshape(3,-1).T
    offsets = np.array([(i,j,k) for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])
    x = ((grid[:,None,:]+offsets[None,:,:])*box.d).reshape(-1,3)
    velocity = np.array([1e5,-5e4,2.5e4])
    v = np.broadcast_to(velocity,x.shape).copy()
    w = np.full(len(x),density*box.cell_volume/8)
    metal = MetalKernels(SimpleNamespace(box=box,lo=np.array(box.lo),d=box.d,n_ref=density))
    order = np.random.default_rng(5).permutation(len(x)); rows=[]
    for label, indices in (('ordered',np.arange(len(x))),('shuffled',order),('partial_group',order[:-3])):
        result = metal.moments(x[indices],v[indices],w[indices])
        expected_number = float(w[indices].sum())
        row = dict(case=label,ions=len(indices),
            number_relative_error=float((result[0].sum()*box.cell_volume-expected_number)/expected_number),
            integrated_flux_relative_error=float(np.max(abs(result[1:].sum(axis=(1,2,3))*box.cell_volume
                -expected_number*velocity))/np.max(abs(expected_number*velocity))))
        if label != 'partial_group':
            row['uniform_density_Linf_relative']=float(np.max(abs(result[0]-density))/density)
        rows.append(row)
    save(out,'cic_analytic',rows);print('analytic CIC',rows,flush=True)


def particle_block(model):
    """Actual device particles, including GPU casts and synchronization."""
    metal = MetalKernels(model)
    cpu = lambda: deposit(model.x, model.v, model.w, model.lo, model.d, model.shape, model.nthreads)
    gpu = lambda: metal.moments(model.x, model.v, model.w)
    reference, measured = cpu(), gpu()
    times = {}
    for name, operation in (('cpu', cpu), ('metal', gpu)):
        samples = []
        for _ in range(3):
            start = time.perf_counter(); operation(); samples.append(time.perf_counter()-start)
        times[name] = float(np.median(samples))
    inputs = metal.prepare_moments(model.x, model.v, model.w)
    metal.mx.eval(*inputs)
    samples = []
    for _ in range(5):
        start = time.perf_counter(); result = metal.moments_device(inputs);metal.mx.eval(result)
        samples.append(time.perf_counter()-start)
    times['metal_resident_kernel_only'] = float(np.median(samples))
    return dict(ions=len(model.x), grid=list(map(int, model.shape)), seconds=times,
        speedup=times['cpu']/times['metal'],
        resident_kernel_only_speedup=times['cpu']/times['metal_resident_kernel_only'],
        density_relative_L2=float(np.linalg.norm(measured[0]-reference[0])/np.linalg.norm(reference[0])),
        density_integral_relative_error=float((measured[0].sum()*model.box.cell_volume-model.w.sum())/model.w.sum()),
        flux_relative_L2=float(np.linalg.norm(measured[1:]-reference[1:])/np.linalg.norm(reference[1:])))


def device(out, dx=.01, dz=.01, steps=40, modes=(0, 1, 2)):
    state = full_flux('outputs/fine_uniform_h2', 0.)
    previous = out/'device.json'
    rows = json.loads(previous.read_text()) if previous.exists() else []
    for gpu in modes:
        start = time.perf_counter()
        model, info = build(state, dx, dz, ppc=8, floor=.05, seed=0,
            drive=2, duration_s=40e-6, dt_ns=.5, thermal=True, gpu=gpu)
        init_s = time.perf_counter()-start
        initial = diagnostics(model)
        run = out/f'gpu_{gpu}'; run.mkdir(exist_ok=True)
        (run/'setup.json').write_text(json.dumps(info, indent=1)+'\n')
        times = []; history = [initial]
        for step in range(steps):
            start = time.perf_counter(); model.step(); model.source.advance(model)
            times.append(time.perf_counter()-start)
            if step == 9 or step == steps-1: history.append(diagnostics(model))
        last = history[-1]
        save(run, 'history', history)
        save_snapshot(model, run, 0, np.random.default_rng(1))
        row = dict(gpu_mode=gpu, initialization_s=init_s, measured_steps=steps-10,
                   measured_step_wall_s=sum(times[10:]), time_us=model.time*1e6,
                   field_substeps=model.substeps, endpoint=last)
        rows.append(row);save(out, 'device', rows)
        print('device',gpu,row['measured_step_wall_s'],last['load_energy_residual_j'],flush=True)
        if not (out/'particle_block.json').exists():
            result = particle_block(model);save(out, 'particle_block', result)
            print('particle block',result,flush=True)
        del model;gc.collect()
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out')
    parser.add_argument('--section', choices=('independent', 'device', 'all'), default='all')
    parser.add_argument('--dx', type=float, default=.01)
    parser.add_argument('--dz', type=float, default=.01)
    parser.add_argument('--steps', type=int, default=40)
    parser.add_argument('--modes', type=int, nargs='+', choices=(0, 1, 2), default=(0, 1, 2))
    args = parser.parse_args();out=Path(args.out);out.mkdir(parents=True, exist_ok=True)
    project = Path(__file__).resolve().parents[1]
    paths = [project/'plasma_bench'/name for name in
        ('hybrid3d_metal.py', 'hybrid3d_metal_energy.py', 'hybrid3d_device.py', 'hybrid3d_energy_electrons.py')]
    paths += [Path(__file__), project/'examples'/'hybrid3d_frc.py']
    save(out, 'source_hashes', {str(p.relative_to(project)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    if args.section in ('independent', 'all'):independent(out)
    if args.section in ('device', 'all'):device(out,args.dx,args.dz,args.steps,args.modes)
