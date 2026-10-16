"""Compare complete biased-coil trajectories by their actual physical response."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as E, electron_mass as ME
from plasma_bench.material_experiment import build, source_configs
from plasma_bench.material_plasma_samples import material_samples
from plasma_bench.material_vacuum_samples import physical_samples
from .material_coil_plot import state_at
from .implicit_plasma_records import sha


def drift_content(response):
    """How much plasma and resistive power occupy large-drift samples?"""
    fields = response['fields']; n = fields['n'].ravel()
    te = fields['temperature'].reshape(-1, 2)[:, 0]
    drift2 = np.sum(response['drift']**2, axis=1)
    ratio = np.sqrt(drift2/(E*te/ME)); volume = response['volume']
    weights = {'volume': volume, 'particle': n*volume,
               'electron_thermal_energy': n*te*volume,
               'current_squared_volume': (E*n)**2*drift2*volume}
    return dict(maximum_drift_over_electron_thermal_speed=float(ratio.max()),
        scope='Actual quadrature; current-squared volume weights equal Joule-power weights only for the prescribed uniform scalar resistivity. Thresholds describe content, not a validity theorem.',
        thresholds=[dict(drift_over_thermal=threshold,
            fractions={name: float(weight[ratio >= threshold].sum()/weight.sum())
                       if weight.sum() else None for name, weight in weights.items()})
            for threshold in (.1, 1.)])


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
    labels = []
    for data in (first, second):
        pulsed = [g for g in source_configs(data['config']) if g['pulse_center_t']['amplitude'] != 0.]
        labels.append('Bias only' if not pulsed else ' + '.join(
            f"{abs(g['pulse_center_t']['amplitude'])*1e3:g} mT · {len(g['coils_radius_z_turns'])} pulse coils" for g in pulsed))
    field_cap = max(abs(x) for path, data in runs for x in data['records'][-1]['plasma_bz_range_t'])*1e6
    increment = 10**(np.floor(np.log10(field_cap))-1) if field_cap else 1.
    field_cap = max(increment, np.ceil(field_cap/increment)*increment)
    for run_index, ((path, data), label, color) in enumerate(zip(runs, labels, colors)):
        rows = [data['initial']]+data['records']; t = np.array([r['time'] for r in rows])*1e6
        f, model, meta = build(data['config'])
        side = f.domain.ids[-1]; midplane_radius = []; rms_radius = []
        for row in rows:
            u = state_at(path, row.get('step', 0))
            midplane_radius.append(float(np.interp(0., model.space.mesh.zf, u.fluid.coordinates[side, 0])))
            points = f.geometry(u.fluid.coordinates)['points']
            rms_radius.append(float(np.sqrt(np.sum(f.particles*points[..., 0]**2)/f.particles.sum())))
        axes[0, 0].plot(t, [r['nominal_vacuum_center_field_t']*1e3 for r in rows], color=color, label=label)
        axes[0, 1].plot(t, [r['maximum_plasma_radius_m']*1e3 for r in rows], color=color, label=label+' · max')
        axes[0, 1].plot(t, np.array(midplane_radius)*1e3, '--', color=color, label=label+' · middle')
        axes[0, 1].plot(t, np.array(rms_radius)*1e3, ':', color=color, label=label+' · particle RMS')
        axes[0, 2].plot(t, [r['joule_heat_j']*1e3 for r in rows], color=color, label=label+' · Joule')
        numerical_heat = [r.get('numerical_viscous_heat_j') for r in rows]
        if all(value is not None for value in numerical_heat):
            axes[0, 2].plot(t, np.array(numerical_heat)*1e3, '--', color=color,
                           label=label+' · numerical ion heat')
        axes[1, 0].plot(t, [r['peak_current_components_a_m2'][1]/1e3 for r in rows], color=color, label=label)
        state = state_at(path, data['completed_steps'])
        r = np.linspace(.001, .999, 120)*data['config']['plasma']['radius_m']
        z = np.linspace(-.999, .999, 360)*data['config']['plasma']['half_length_m']
        rr, zz = np.meshgrid(r, z, indexing='ij')
        n0 = meta['central_density_m3']*((1-(rr/data['config']['plasma']['radius_m'])**2)*(1-(zz/data['config']['plasma']['half_length_m'])**2))**1.5
        samples = material_samples(model, state, rr, zz, n0, time=data['duration'])
        axis = axes[1, run_index+1]
        artist = axis.pcolormesh(samples['points'][..., 1]*100, samples['points'][..., 0]*100,
            samples['magnetic'][..., 2]*1e6, shading='auto', cmap='RdBu_r', vmin=-field_cap, vmax=field_cap, rasterized=True)
        axis.contour(samples['points'][..., 1]*100, samples['points'][..., 0]*100,
            samples['flux'], levels=12, colors='#203c48', linewidths=.55)
        axis.set(xlabel='Axial position z [cm]', ylabel='Radius r [cm]', title=f"{label}: final Bz and flux\n{data['duration']*1e6:g} µs; radial scale expanded")
        fig.colorbar(artist, ax=axis, label='Axial field Bz [µT]', shrink=.8)
        record['runs'][path.name] = dict(source_sha256=data['source_sha256'], summary_sha256=sha(path/'summary.json'),
            label=label, particle_rms_radius_history_m=rms_radius,
            midplane_radius_history_m=midplane_radius, time_history_s=[r['time'] for r in rows],
            joule_heat_history_j=[r['joule_heat_j'] for r in rows],
            numerical_viscous_heat_history_j=numerical_heat,
            numerical_bulk_viscosity=data['config'].get('numerical_bulk_viscosity'),
            maximum_absolute_external_work_residual_over_initial_energy=max(
                abs(r['relative_external_work_residual']) for r in rows),
            minimum_midplane_radius_m=min(midplane_radius),
            minimum_particle_rms_radius_m=min(rms_radius),
            final_midplane_radius_m=midplane_radius[-1],
            final={key: rows[-1][key] for key in ('time', 'maximum_plasma_radius_m', 'maximum_plasma_half_length_m',
                'plasma_bz_range_t', 'peak_current_components_a_m2', 'joule_heat_j', 'external_work_j',
                'particle_mean_temperature_ev', 'maximum_drift_over_electron_thermal_speed')})
        record['runs'][path.name]['final']['numerical_viscous_heat_j'] = numerical_heat[-1]
        # Compare independently reconstructed exterior coil flux with the
        # volume-current linkage, using the actual peak-drive saved state.
        peak = max(rows, key=lambda row: abs(row['nominal_vacuum_center_field_t']))
        peak_state = state_at(path, peak.get('step', 0))
        response = model.response(peak_state, rate=False, time=peak['time'])
        record['runs'][path.name]['drift_content'] = {
            'peak_drive': dict(time_s=peak['time'], **drift_content(response)),
            'final': dict(time_s=data['duration'], **drift_content(
                model.response(state, rate=False, time=data['duration'])))}
        group_comparisons = []
        for index, group in enumerate(meta['source_groups']):
            coils = np.asarray(group['coils_radius_z_turns'])
            flux, _, _ = physical_samples(model.vacuum, peak_state.fluid.coordinates,
                model.space.vacuum_trace@peak_state.magnetic, coils[:, 0], coils[:, 1],
                geometry=response['vacuum_response'])
            exterior = float(2*np.pi*np.sum(coils[:, 2]*flux))
            volume = float(response['source_group_linkages_wb'][index])
            group_comparisons.append(dict(name=group['name'], time_s=peak['time'],
                exterior_coil_flux_linkage_wb=exterior, volume_current_linkage_wb=volume,
                relative_difference=exterior/volume-1 if volume else None))
        record['runs'][path.name]['peak_drive_linkage_comparison'] = group_comparisons
    labels = [('Prescribed field at centre', 'Vacuum centre Bz [mT]'),
              ('Compression at the middle versus outer expansion', 'Plasma radius [mm]'),
              ('Resistive and numerical compression heating', 'Accumulated heat [mJ]'), ('Azimuthal plasma current', 'Peak |Jθ| [kA/m²]')]
    for axis, (title, unit) in zip([*axes[0], axes[1, 0]], labels):
        axis.set(title=title, xlabel='Time [µs]', ylabel=unit); axis.grid(alpha=.2); axis.legend(frameon=False, fontsize=8)
    a, b = [d['records'][-1] for p, d in runs]
    record['final_radius_second_minus_first_m'] = b['maximum_plasma_radius_m']-a['maximum_plasma_radius_m']
    tube=second['config']['tube'];plasma=second['config']['plasma']
    fig.suptitle(f"Coil-pulse response: magnetic field, current, heat and plasma motion\n{tube['length_m']:g} m × {tube['radius_m']*200:g} cm tube · assumed prepared {plasma['gas']} at {plasma['neutral_fill_mtorr']:g} mTorr inventory · no neutral/ionization model", fontsize=14)
    fig.savefig(out/'pulse_vs_bias.png', dpi=170); fig.savefig(out/'pulse_vs_bias.pdf'); plt.close(fig)
    (out/'comparison.json').write_text(json.dumps(record, indent=2)+'\n')
    print('Physical coil comparison complete:', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); parser.add_argument('control'); parser.add_argument('pulse')
    args = parser.parse_args(); compare(args.destination, args.control, args.pulse)
