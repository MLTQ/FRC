"""Passive CIC power/position-moment comparison on saved native particles."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
from types import SimpleNamespace

if __name__ == '__main__':
    os.environ.setdefault('NUMBA_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np
from numba import get_num_threads
from scipy.constants import elementary_charge as qe
from plasma_bench.hybrid3d import deposit
from plasma_bench.hybrid3d_device import Device, device_box
from plasma_bench.hybrid3d_metal import MetalKernels
from plasma_bench.hybrid3d_metal_cic import SplitCoordinateMetal
import examples.hybrid3d_gpu_check as analytic_check


def run(args):
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    hashes = {}
    for relative in ('examples/hybrid3d_cic_work.py', 'plasma_bench/hybrid3d_metal_cic.py',
                     'plasma_bench/hybrid3d_metal.py', 'plasma_bench/hybrid3d.py'):
        for path in (root/relative, (root/relative).with_suffix('.md')):
            copy = out/'source'/path.relative_to(root)
            copy.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, copy)
            hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    info = json.loads(Path(args.setup).read_text())
    d = info['device']
    device = Device(wall_radius=d['wall'], coil_radius=d['coil'], tube_half_length=d['half'],
                    end_gap=d['end_gap'], inner_radius=d['inner'])
    box = device_box(device, args.dx, args.dz)
    with np.load(args.state) as state:
        x, v, w = state['x'], state['v'], state['w']
        if list(box.shape) != list(state['n'].shape):
            raise ValueError('Reconstructed geometry does not match native snapshot')
        time_s = float(state['time_s'])
    cells = np.floor((x-box.lo)/box.d)
    if np.any(cells < 0) or np.any(cells+1 >= np.asarray(box.shape)):
        raise ValueError('Linear-moment comparison requires no periodic stencil crossing')
    # Electrostatic vacuum quadrupole E=(-x,-y,2z) V/m, coefficient 1 V/m².
    # Trilinear interpolation reproduces this exact linear field.
    analytic_e = x*np.array([-1., -1., 2.])
    individual = qe*w*np.sum(v*analytic_e, axis=1)
    exact_power = float(individual.sum())
    throughput = float(np.abs(individual).sum())
    expected_number = float(w.sum())
    expected_position_moment = np.sum(w[:, None]*x, axis=0)
    nodes = box.nodes(0)[:, None, None], box.nodes(1)[None, :, None], box.nodes(2)[None, None, :]
    model = SimpleNamespace(box=box, lo=np.asarray(box.lo), d=box.d, n_ref=info['n_peak'])
    cpu = lambda: deposit(x, v, w, model.lo, model.d, np.asarray(box.shape), get_num_threads())
    reference = cpu()
    rows = []
    for name, kernel in (('cpu64', None), ('world_coordinate_metal', MetalKernels(model)),
                         ('split_coordinate_metal', SplitCoordinateMetal(model))):
        operation = cpu if kernel is None else lambda: kernel.moments(x, v, w)
        operation()  # compile/warm first; timing below includes host handling.
        times, repeats = [], []
        first = None
        for _ in range(3):
            start = time.perf_counter()
            moments = operation()
            times.append(time.perf_counter()-start)
            if first is None:
                first = moments
            else:
                repeats.append(float(np.linalg.norm(moments-first)/np.linalg.norm(first)))
        power = qe*box.cell_volume*sum(s*np.sum(first[a+1]*nodes[a]) for a,s in enumerate((-1.,-1.,2.)))
        moment = np.array([np.sum(first[0]*node)*box.cell_volume for node in nodes])
        rows.append(dict(kernel=name, median_seconds=float(np.median(times)),
            density_relative_l2_to_cpu=float(np.linalg.norm(first[0]-reference[0])/np.linalg.norm(reference[0])),
            flux_relative_l2_to_cpu=float(np.linalg.norm(first[1:]-reference[1:])/np.linalg.norm(reference[1:])),
            number_relative_error=float((first[0].sum()*box.cell_volume-expected_number)/expected_number),
            position_moment_error_m_per_ion=((moment-expected_position_moment)/expected_number).tolist(),
            quadrupole_power_w=float(power), exact_power_w=exact_power,
            power_absolute_error_w=float(abs(power-exact_power)),
            power_relative_error_to_net=float(abs(power-exact_power)/abs(exact_power)) if exact_power else None,
            power_relative_error_to_sum_abs=float(abs(power-exact_power)/throughput),
            repeat_relative_l2=repeats))
        report = dict(time_s=time_s, ions=len(x), rows=rows, source_hashes=hashes,
                      native_state=str(Path(args.state).resolve()),
                      state_sha256=hashlib.sha256(Path(args.state).read_bytes()).hexdigest(),
                      geometry=dict(shape=list(box.shape), lo=list(box.lo), hi=list(box.hi)),
                      equations='E=(-x,-y,2z), P=sum(q w v dot E(x))=integral(J dot E)dV; no periodic crossing',
                      note='Passive manufactured-field evaluation on actual saved particles; this field is never inserted into a physical run.')
        (out/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False))
        print(rows[-1], flush=True)
    for label, cls in (('original',MetalKernels),('split',SplitCoordinateMetal)):
        directory = out/label
        directory.mkdir()
        analytic_check.MetalKernels = cls
        analytic_check.independent_cic(directory)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out')
    parser.add_argument('--state', required=True)
    parser.add_argument('--setup', required=True)
    parser.add_argument('--dx', type=float, default=.01)
    parser.add_argument('--dz', type=float, default=.01)
    run(parser.parse_args())
