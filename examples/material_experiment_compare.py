"""Compare complete biased-coil trajectories by their actual physical response."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from plasma_bench.material_experiment import build
from plasma_bench.material_plasma_samples import material_samples
from .material_coil_plot import state_at
from .implicit_plasma_records import sha


def compare(destination, control_path, pulse_path):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh physical comparison directory required')
    runs = []
    for source in (control_path, pulse_path):
        path = Path(source); data = json.loads((path/'summary.json').read_text())
        if not data['terminal'] or not data['success']:
            raise ValueError('Complete successful trajectories required')
        runs.append((path, data))
    first, second = [data for path, data in runs]
    for key in ('tube', 'plasma', 'mesh', 'evolution'):
        if first['config'][key] != second['config'][key]:
            raise ValueError('Matching chamber, plasma and evolution required')
    out.mkdir(); fig, axes = plt.subplots(2, 3, figsize=(15, 8.6), constrained_layout=True)
    record = {'scope': 'Actual saved trajectories; exploratory pre-ionized plasma, not converged formation or confinement.', 'runs': {}}
    colors = ['#4b7281', '#b84b33']
    amplitude = abs(second['config']['source']['pulse_center_t']['amplitude'])*1e3
    labels = ('Bias only', f'{amplitude:g} mT opposing pulse')
    field_cap = max(abs(x) for path, data in runs for x in data['records'][-1]['plasma_bz_range_t'])*1e6
    increment = 10**(np.floor(np.log10(field_cap))-1) if field_cap else 1.
    field_cap = max(increment, np.ceil(field_cap/increment)*increment)
    for (path, data), label, color in zip(runs, labels, colors):
        rows = [data['initial']]+data['records']; t = np.array([r['time'] for r in rows])*1e6
        f, model, meta = build(data['config'])
        side = f.domain.ids[-1]; midplane_radius = []
        for row in rows:
            u = state_at(path, row.get('step', 0))
            midplane_radius.append(float(np.interp(0., model.space.mesh.zf, u.fluid.coordinates[side, 0])))
        axes[0, 0].plot(t, [r['nominal_vacuum_center_field_t']*1e3 for r in rows], color=color, label=label)
        axes[0, 1].plot(t, [r['maximum_plasma_radius_m']*1e3 for r in rows], color=color, label=label+' · max')
        axes[0, 1].plot(t, np.array(midplane_radius)*1e3, '--', color=color, label=label+' · middle')
        axes[0, 2].plot(t, [r['joule_heat_j']*1e6 for r in rows], color=color, label=label)
        axes[1, 0].plot(t, [r['peak_current_components_a_m2'][1]/1e3 for r in rows], color=color, label=label)
        state = state_at(path, data['completed_steps'])
        r = np.linspace(.001, .999, 120)*data['config']['plasma']['radius_m']
        z = np.linspace(-.999, .999, 360)*data['config']['plasma']['half_length_m']
        rr, zz = np.meshgrid(r, z, indexing='ij')
        n0 = meta['central_density_m3']*((1-(rr/data['config']['plasma']['radius_m'])**2)*(1-(zz/data['config']['plasma']['half_length_m'])**2))**1.5
        samples = material_samples(model, state, rr, zz, n0, time=data['duration'])
        axis = axes[1, 1] if label == labels[0] else axes[1, 2]
        artist = axis.pcolormesh(samples['points'][..., 1]*100, samples['points'][..., 0]*100,
            samples['magnetic'][..., 2]*1e6, shading='auto', cmap='RdBu_r', vmin=-field_cap, vmax=field_cap, rasterized=True)
        axis.contour(samples['points'][..., 1]*100, samples['points'][..., 0]*100,
            samples['flux'], levels=12, colors='#203c48', linewidths=.55)
        axis.set(xlabel='Axial position z [cm]', ylabel='Radius r [cm]', title=f"{label}: final Bz and flux\n{data['duration']*1e6:g} µs; radial scale expanded")
        fig.colorbar(artist, ax=axis, label='Axial field Bz [µT]', shrink=.8)
        record['runs'][path.name] = dict(source_sha256=data['source_sha256'], summary_sha256=sha(path/'summary.json'),
            midplane_radius_history_m=midplane_radius, time_history_s=[r['time'] for r in rows],
            final_midplane_radius_m=midplane_radius[-1],
            final={key: rows[-1][key] for key in ('time', 'maximum_plasma_radius_m', 'maximum_plasma_half_length_m',
                'plasma_bz_range_t', 'peak_current_components_a_m2', 'joule_heat_j', 'external_work_j',
                'particle_mean_temperature_ev', 'maximum_drift_over_electron_thermal_speed')})
    labels = [('Prescribed field at centre', 'Vacuum centre Bz [mT]'),
              ('Compression at the middle versus outer expansion', 'Plasma radius [mm]'),
              ('Accumulated Joule heat', 'Heat [µJ]'), ('Azimuthal plasma current', 'Peak |Jθ| [kA/m²]')]
    for axis, (title, unit) in zip([*axes[0], axes[1, 0]], labels):
        axis.set(title=title, xlabel='Time [µs]', ylabel=unit); axis.grid(alpha=.2); axis.legend(frameon=False, fontsize=8)
    a, b = [d['records'][-1] for p, d in runs]
    record['final_radius_pulse_minus_control_m'] = b['maximum_plasma_radius_m']-a['maximum_plasma_radius_m']
    tube=second['config']['tube'];plasma=second['config']['plasma']
    fig.suptitle(f"Coil-pulse response: magnetic field, current, heat and plasma motion\n{tube['length_m']:g} m × {tube['radius_m']*200:g} cm tube · assumed prepared {plasma['gas']} at {plasma['neutral_fill_mtorr']:g} mTorr inventory · no neutral/ionization model", fontsize=14)
    fig.savefig(out/'pulse_vs_bias.png', dpi=170); fig.savefig(out/'pulse_vs_bias.pdf'); plt.close(fig)
    (out/'comparison.json').write_text(json.dumps(record, indent=2)+'\n')
    print('Physical pulse/control comparison complete:', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); parser.add_argument('control'); parser.add_argument('pulse')
    args = parser.parse_args(); compare(args.destination, args.control, args.pulse)
