"""Strict nonlinear toroidal record audit, comparisons and scientific figures."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
from .toroidal_cases import protocol
from .toroidal_measurements import provenance
from .verification_records import write_record


def validate(directory):
    data = json.loads((directory/'measurements.json').read_text())
    cases = protocol()
    if not data['completed'] or data['provenance'] != provenance() or [r['name'] for r in data['cases']] != [c.name for c in cases]:
        raise ValueError('Incomplete or stale toroidal protocol')
    for row, case in zip(data['cases'], cases):
        if not row['completed'] or row['configuration'] != case.configuration() or row['provenance'] != data['provenance']:
            raise ValueError('Toroidal case/source identity differs')
        if json.loads((directory/(case.name+'.json')).read_text()) != row:
            raise ValueError('Collected and individual records differ')
        path = directory/row['fields']['file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['fields']['sha256']:
            raise ValueError('Lossless toroidal archive hash differs')
        with np.load(path, allow_pickle=False) as fields:
            if any(not np.isfinite(fields[key]).all() for key in fields.files):
                raise ValueError('Nonfinite archived toroidal fields')
        history = row['history_columns']
        if any(len(values) != row['final']['steps']+1 for values in history.values()):
            raise ValueError('Incomplete accepted-step history')
        if min(history['physical_ohmic_heat_j']) < 0 or min(history['min_pressure_pa']) <= 0 or min(history['min_density_kg_m3']) <= 0:
            raise ValueError('Nonphysical accepted diagnostics')
        if abs(history['time_s'][-1]-case.end_time) > 1e-13*case.end_time:
            raise ValueError('Case did not reach declared final time')
    return data


def time_comparisons(directory, data):
    rows = [r for r in data['cases'] if r['configuration']['control'].startswith('time')]
    rows.sort(key=lambda r: r['configuration']['dt_max'], reverse=True)
    with np.load(directory/rows[-1]['fields']['file'], allow_pickle=False) as field:
        reference, volume = field['final_conserved'].copy(), field['volumes'].copy()
    errors = []
    for row in rows[:-1]:
        with np.load(directory/row['fields']['file'], allow_pickle=False) as field:
            difference = field['final_conserved']-reference
        # Scale each conserved variable by its reference RMS to keep dimensions separate.
        values = {}
        for name, index in (('density', 0), ('axial_momentum', 2), ('energy', 3), ('angular_momentum', 4)):
            values[name] = float(np.sqrt(np.sum(volume*difference[..., index]**2)/np.sum(volume*reference[..., index]**2)))
        errors.append({'name': row['name'], 'dt_max_s': row['configuration']['dt_max'], 'relative_volume_errors': values})
    ratios = {key: [errors[i]['relative_volume_errors'][key]/errors[i+1]['relative_volume_errors'][key] for i in (0, 1)] for key in errors[0]['relative_volume_errors']}
    return {'reference_case': rows[-1]['name'], 'comparisons': errors, 'halving_ratios': ratios}


def plots(directory, data, times):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    fig.suptitle('Nonlinear toroidal MHD: measured transport, waves, resistive decay and conservation', fontsize=14)
    for kind, axis in (('contact', axes[0, 0]), ('wave', axes[0, 1])):
        rows = [r for r in data['cases'] if r['configuration']['kind'] == kind and r['configuration']['control'] == 'base']
        for index, (key, label) in enumerate(((('density_relative_L2' if kind == 'contact' else 'rotation_relative_L2'), 'Density' if kind == 'contact' else 'Rotation'),
                           (('rotation_relative_L2' if kind == 'contact' else 'toroidal_field_relative_L2'), 'Rotation' if kind == 'contact' else 'Toroidal field'))):
            physical_key = ('density' if index == 0 else 'rotation_velocity') if kind == 'contact' else ('rotation_velocity' if index == 0 else 'toroidal_field')
            axis.loglog([r['configuration']['nz'] for r in rows], [r['errors'][key] for r in rows], 'o-', color=f'C{index}', label=label+' means')
            axis.loglog([r['configuration']['nz'] for r in rows], [r['errors']['physical_reconstruction_relative_L2'][physical_key] for r in rows], 's--', color=f'C{index}', label=label+' reconstructed')
        axis.set(title='Exact finite-amplitude rotating transport' if kind == 'contact' else 'Linear torsional Alfven wave', xlabel='Axial cells', ylabel='Relative volume L2 error')
        axis.legend(fontsize=8)
        axis.set_xticks([24, 48, 96, 192])
    rows = [r for r in data['cases'] if r['configuration']['kind'] == 'decay']
    for key, label in (('toroidal_field_relative_L2', 'Field cell mean'), ('axis_g_relative_error', 'Axis g'), ('physical_heat_relative_error', 'Physical heat')):
        axes[0, 2].loglog([r['configuration']['nr'] for r in rows], [abs(r['errors'][key]) for r in rows], 'o-', label=label)
    axes[0, 2].loglog([r['configuration']['nr'] for r in rows], [r['errors']['physical_reconstruction_relative_L2']['toroidal_field'] for r in rows], 's--', label='Reconstructed field')
    axes[0, 2].set(title='Independent conducting-Bessel decay', xlabel='Radial cells', ylabel='Absolute relative error')
    axes[0, 2].set_xticks([8, 16, 32, 64, 128])
    axes[0, 2].legend(fontsize=9)
    for key, label in (('density', 'Density'), ('angular_momentum', 'Angular momentum')):
        axes[1, 0].loglog([r['dt_max_s'] for r in times['comparisons']], [r['relative_volume_errors'][key] for r in times['comparisons']], 'o-', label=label)
    axes[1, 0].set(title='Time controls against finer same-grid run', xlabel='Maximum step (s)', ylabel='Relative difference')
    axes[1, 0].legend(fontsize=9)
    row = next(r for r in data['cases'] if r['name'] == 'contact_nr8_nz192_base')
    for key, label in (('mass_kg', 'Mass'), ('angular_momentum_kg_m2_s', 'Angular momentum'), ('total_energy_j', 'Total energy')):
        values = np.asarray(row['history_columns'][key])
        axes[1, 1].plot(row['history_columns']['time_s'], 1e15*(values-values[0])/values[0], label=label)
    axes[1, 1].set(title='Actual nonlinear conserved budgets', xlabel='Time (s)', ylabel='Relative change (×10⁻¹⁵)')
    axes[1, 1].legend(fontsize=9)
    for tag in ('rotation', 'pinch', 'combined'):
        rows = [r for r in data['cases'] if r['configuration']['kind'] == 'equilibrium' and r['configuration']['control'] == tag]
        axes[1, 2].semilogy([r['configuration']['nr'] for r in rows], [r['errors']['radial_velocity_rms_m_s'] for r in rows], 'o-', label=tag)
    axes[1, 2].set(title='Exact equilibrium spurious radial flow', xlabel='Radial cells', ylabel='Volume RMS velocity (m/s)')
    axes[1, 2].legend(fontsize=9)
    axes[1, 2].set_xticks([8, 16, 32])
    for ax in axes.flat:
        ax.grid(alpha=.2)
    for ax in axes[0]:
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    fig.savefig(directory/'toroidal_accuracy.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
    fig.suptitle('Finite-amplitude rotating entropy wave: exact versus evolved at t = 0.4 s')
    for nz in (24, 192):
        row = next(r for r in data['cases'] if r['name'] == f'contact_nr8_nz{nz}_base')
        with np.load(directory/row['fields']['file'], allow_pickle=False) as field:
            z = .5*(field['z_faces'][1:]+field['z_faces'][:-1])
            for ax, index, label in ((axes[0], 0, 'Density (kg/m³)'), (axes[1], 4, 'Angular momentum density (kg/m/s)')):
                ax.plot(z, field['final_conserved'][4, :, index], '.-', label=f'{nz} axial cells', markersize=3)
                if nz == 192:
                    ax.plot(z, field['reference_conserved'][4, :, index], 'k--', label='Exact cell averages', linewidth=1)
                ax.set(xlabel='z (m)', ylabel=label)
                ax.grid(alpha=.2)
    for ax in axes:
        ax.legend(fontsize=9)
    fig.savefig(directory/'rotating_transport.png', dpi=160)
    plt.close(fig)


def run(directory='outputs/toroidal_verification'):
    directory = Path(directory)
    data = validate(directory)
    times = time_comparisons(directory, data)
    plots(directory, data, times)
    metrics = {'completed': True, 'case_count': len(data['cases']), 'source_sha256': data['provenance']['source_sha256'], 'time_controls': times,
               'maximum_relative_mass_change': max(abs(r['final']['mass_kg_change'])/r['initial']['mass_kg'] for r in data['cases']),
               'maximum_relative_total_energy_change': max(abs(r['final']['total_energy_j_change'])/r['initial']['total_energy_j'] for r in data['cases']),
               'maximum_physical_error_quadrature_change': max(r['errors']['physical_error_quadrature_absolute_change'] for r in data['cases']),
               'maximum_equilibrium_radial_rms_m_s': max(r['errors']['radial_velocity_rms_m_s'] for r in data['cases'] if r['configuration']['kind'] == 'equilibrium'),
               'report_source_files': {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in ('examples/verify_toroidal.py', 'examples/verification_records.py')},
               'figures': {name: hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in ('toroidal_accuracy.png', 'rotating_transport.png')}}
    write_record(directory/'metrics.json', metrics)
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    run()
