"""Dynamic vacuum-mask startup replay (frc-ibl.11.9; see vacuum_mask_replay.md).

    python -m examples.vacuum_mask_replay run OUT seed=0 mask=1 min_ions=1 t_end_us=1 [ppc=8]
    python -m examples.vacuum_mask_replay compare BASE0 BASE1 MASK0 MASK1

`run` repeats the v4 baseline shot (outputs/fine_uniform_h2 from t = 0; 20 x 20 x 40 mm cells, 8 ions/cell, 5% floor,
coupled banks, energy electrons, 1 ns) through examples.hybrid3d_frc.main unchanged, with the model class replaced
by VacuumMaskHybrid when mask=1 (mask=0 is the plain baseline, for a second seed). Mask counters are added to
every history row. `ppc` changes only the macro-ions per cell (the ion-count convergence replay, mask=0).

`compare` tabulates the counters of the four runs (two seeds, without and with the mask) and, at shared snapshot
times, the seed-to-seed difference of B_response and n on plasma (mean n >= floor) and thin (0 < mean n < floor)
gas nodes.
"""
import functools
import json
import sys
from pathlib import Path
import numpy as np
import examples.hybrid3d_frc as frc
from plasma_bench.hybrid3d import faces_to_nodes
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, GAS
from plasma_bench.hybrid3d_vacuum_mask import VacuumMaskHybrid

SOURCE = 'outputs/fine_uniform_h2'
BASELINE = dict(t_us=0., dx=.02, dz=.04, ppc=8, floor=.05, every_us=.25, snap_us=.1, drive=2, thermal=1, dt_ns=1.)


def run(out, seed=0, mask=1, min_ions=1, t_end_us=1., ppc=8):
    if mask:
        frc.DeviceHybrid = functools.partial(VacuumMaskHybrid, min_ions=int(min_ions))
        plain = frc.diagnostics

        def diagnostics(model):
            row = plain(model)
            row.update(model.mask_diagnostics())
            return row
        frc.diagnostics = diagnostics
    history = frc.main(out, SOURCE, t_end_us=float(t_end_us), seed=int(seed), **dict(BASELINE, ppc=int(ppc)))
    setup = Path(out)/'setup.json'
    info = json.loads(setup.read_text())
    info['vacuum_mask'] = dict(enabled=bool(mask), min_ions=int(min_ions),
                               rule='gas cells with < min_ions macro-ions at both ends of an ion step are field vacuum')
    setup.write_text(json.dumps(info, indent=1))
    return history


COUNTERS = ['magnetic_j', 'response_magnetic_j', 'kinetic_j', 'electron_internal_j', 'ohmic_heat_j', 'cutoff_removed_heat_j',
            'electron_empty_node_energy_j', 'electron_boundary_enthalpy_j', 'mask_projection_response_j', 'load_energy_work_j',
            'load_energy_residual_j', 'trapped_flux_wb', 'separatrix_mid_m', 'closed', 'masked_gas_cells', 'substeps']


def node_classes(run, n_mean, floor):
    info = json.loads((Path(run)/'setup.json').read_text())
    dv = info['device']
    dev = Device(dv['wall'], dv['half'], dv['coil'], dv['end_gap'], dv.get('inner', 0.))
    f = (cell_types(device_box(dev, info['dx'], info['dz']), dev) == GAS).astype(float)
    for a in range(3):
        f = .5*(f+np.roll(f, 1, axis=a))
    live = f > 0
    return dict(plasma=live & (n_mean >= floor), thin=live & (n_mean > 0) & (n_mean < floor), empty=live & (n_mean <= 0))


def seed_difference(a, b, floor):
    """Relative seed-to-seed difference of node B_response and n per node class, at shared snapshot times."""
    rows = []
    for path in sorted(Path(a).glob('snap_*.npz')):
        other = Path(b)/path.name
        if not other.exists():
            continue
        sa, sb = np.load(path), np.load(other)
        if abs(float(sa['time_s'])-float(sb['time_s'])) > 1e-10:      # final snapshots land at different times
            continue
        na, nb = sa['n'].astype(float), sb['n'].astype(float)
        ba = np.asarray(faces_to_nodes(*sa['b_response'].astype(float)))
        bb = np.asarray(faces_to_nodes(*sb['b_response'].astype(float)))
        row = dict(time_us=round(float(sa['time_s'])*1e6, 3))
        for name, m in node_classes(a, .5*(na+nb), floor).items():
            if not m.any():
                continue
            db = np.sqrt(np.sum((ba-bb)[:, m]**2)/max(np.sum((.5*(ba+bb))[:, m]**2), 1e-300))
            dn = np.sum(np.abs(na-nb)[m])/max(np.sum((.5*(na+nb))[m]), 1e-300)
            row[name] = dict(nodes=int(m.sum()), b_rel_rms=float(db), n_rel_l1=float(dn))
        rows.append(row)
    return rows


def compare(base0, base1, mask0, mask1):
    runs = dict(base0=base0, base1=base1, mask0=mask0, mask1=mask1)
    table = {}
    end = json.loads((Path(mask0)/'history.json').read_text())[-1]['time_s']
    for name, run in runs.items():
        h = json.loads((Path(run)/'history.json').read_text())
        table[name] = [dict(time_us=round(r['time_s']*1e6, 3), **{k: r.get(k) for k in COUNTERS})
                       for r in h if r['time_s'] <= end+1e-12 and abs(r['time_s']*4e6-round(r['time_s']*4e6)) < 1e-6]
        # the 40 us baseline is cut to the replay window; rows on the 0.25 us marks only
    floor = json.loads((Path(base0)/'setup.json').read_text())['floor_m3']
    out = dict(runs=runs, counters=table,
               seed_difference=dict(baseline=seed_difference(base0, base1, floor), mask=seed_difference(mask0, mask1, floor)))
    Path(mask0, 'vacuum_mask_compare.json').write_text(json.dumps(out, indent=1))
    for name, rows in table.items():
        r = rows[-1]
        print(name, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
    for name, rows in out['seed_difference'].items():
        print(name, json.dumps(rows[-1]))
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        args = dict(a.split('=', 1) for a in sys.argv[3:])
        run(sys.argv[2], **{k: float(v) for k, v in args.items()})
    else:
        compare(*sys.argv[2:6])
