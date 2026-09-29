"""Audit and render the moving adiabatic-species verification evidence."""
from pathlib import Path
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .warm_audit import audit
from .warm_controls import protocol as control_protocol, source_hash as control_source_hash
from .verification_records import write_record


def time_comparisons(records, arrays):
    results = {}
    for kind in ('contact', 'simple_wave'):
        fine = arrays[f'{kind}_time_2048']['final_electron_pressure']
        base = records[f'{kind}_time_2048']['configuration']['reference']
        errors = []
        for steps in (128, 256, 512):
            values = arrays[f'{kind}_time_{steps}']['final_electron_pressure']
            errors.append(float(np.sqrt(np.mean((values-fine)**2))))
        results[kind] = dict(pressure_rms_difference_pa=errors, ratios=[errors[0]/errors[1], errors[1]/errors[2]],
                             reference_steps=2048, periods=.05, relative_to_uniform_pressure=[e/(base['number_density']*base['electron_ev']*1.602176634e-19) for e in errors])
    return results


def figure_convergence(directory, records, audits):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    resolutions = (16, 32, 64, 128, 256)
    for ax, kind, title in zip(axes, ('contact', 'simple_wave'), ('One full translation', 'Nonlinear compression: half transit')):
        for key, label, color in (('electron_pressure', 'Electron pressure, Q1', '#a144c7'), ('reconstructed_ion_pressure', 'Ion pressure, reconstructed', '#087e8b')):
            errors = [audits[f'{kind}_{nz}']['physical_norms'][key]['perturbation_relative_L2'] for nz in resolutions]
            ax.loglog(resolutions, errors, 'o-', label=label, color=color)
        node = [records[f'{kind}_{nz}']['metrics']['electron_node_perturbation_relative_L2'] for nz in resolutions]
        ax.loglog(resolutions, node, '--', color='#a144c7', alpha=.6, label='Electron node samples')
        anchor = audits[f'{kind}_32']['physical_norms']['electron_pressure']['perturbation_relative_L2']
        ax.loglog(resolutions, anchor*(32/np.array(resolutions))**2, ':', color='gray', label='Second-order guide')
        ax.set(title=title, xlabel='Axial cells', ylabel='Relative L2 error in pressure variation')
        ax.set_xticks(resolutions, labels=[str(n) for n in resolutions])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(which='both', alpha=.15); ax.legend(fontsize=8)
    fig.suptitle('Moving electron/ion verification — ideal common velocity, no Hall or collisions')
    path = directory/'warm_convergence.png'; fig.savefig(path, dpi=180); plt.close(fig)
    return path


def figure_profiles(directory, records, arrays):
    from .warm_reference import WarmReference
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.7), constrained_layout=True)
    for column, kind in enumerate(('contact', 'simple_wave')):
        name = f'{kind}_64'; record, fields = records[name], arrays[name]
        ref = WarmReference(**record['configuration']['reference'])
        z = fields['axial_faces']; centers = .5*(z[1:]+z[:-1])
        fine_z = np.linspace(z[0], z[-1], 800)
        target = ref.fields(fine_z, record['actual_time_s'])
        axes[0, column].plot(fine_z, target['electron_pressure'], color='#a144c7', label='Exact electron')
        axes[0, column].plot(z[:-1], fields['final_electron_pressure'][2], '.', color='#a144c7', label='Electron nodes')
        axes[1, column].plot(fine_z, target['ion_pressure'], color='#087e8b', label='Exact ion')
        axes[1, column].plot(centers, fields['final_ion_pressure'][2], '.', color='#087e8b', label='Ion cell means')
        axes[0, column].set_title('Translation' if kind == 'contact' else 'Nonlinear compression')
        for row in range(2):
            axes[row, column].set(xlabel='z [m]', ylabel='Pressure [Pa]')
            axes[row, column].grid(alpha=.15); axes[row, column].legend(fontsize=8)
    fig.suptitle('Separate species at nz=64 — total energy evolves through conservative fluid flux')
    path = directory/'warm_species_profiles.png'; fig.savefig(path, dpi=180); plt.close(fig)
    return path


def run(directory='outputs/warm_adiabatic_verification'):
    directory = Path(directory)
    data, records, arrays, audits = audit(directory)
    extra_data, extra_records, extra_arrays, extra_audits = audit('outputs/warm_adiabatic_controls', control_protocol())
    assert extra_data['control_source_sha256'] == control_source_hash()
    assert all(r['control_source_sha256'] == control_source_hash() for r in extra_records.values())
    records.update(extra_records); arrays.update(extra_arrays); audits.update(extra_audits)
    times = time_comparisons(records, arrays)
    controls = {}
    for kind in ('contact', 'simple_wave'):
        a, b = (arrays[f'{kind}{suffix}']['final_electron_pressure'] for suffix in ('_128', '_half_cfl_128'))
        reference = records[f'{kind}_128']['configuration']['reference']
        rms = float(np.sqrt(np.mean((a-b)**2)))
        controls[kind+'_half_cfl'] = dict(electron_pressure_rms_difference_pa=rms, relative_to_base_node_error=rms/records[f'{kind}_128']['metrics']['electron_node_error_pa'])
    a, b = arrays['contact_64']['final_electron_pressure'], arrays['contact_scale_64']['final_electron_pressure']
    controls['density_scale'] = dict(relative_pressure_difference=float(np.linalg.norm(a-b)/np.linalg.norm(a)))
    a, b = arrays['contact_64']['final_electron_pressure'][0], arrays['contact_radial_8']['final_electron_pressure'][0]
    controls['radial_resolution'] = dict(relative_pressure_difference=float(np.linalg.norm(a-b)/np.linalg.norm(a)))
    a, b = (arrays[f'contact_common_time_nr{nr}']['final_electron_pressure'][0] for nr in (4, 8))
    controls['radial_resolution_common_time'] = dict(relative_pressure_difference=float(np.linalg.norm(a-b)/np.linalg.norm(a)))
    figures = [figure_convergence(directory, records, audits), figure_profiles(directory, records, arrays)]
    sources = ['examples/'+name+'.py' for name in ('verify_warm_adiabatic', 'warm_audit', 'warm_norms', 'compatible_hall_audit', 'warm_controls')]
    report = dict(completed=True, provenance=data['provenance'], cases=audits, time_comparisons=times, controls=controls,
                  report_sources={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources},
                  figures={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in figures},
                  scope='Moving equal-gamma ideal species verification; not coupled warm Hall or FRC formation')
    write_record(directory/'report.json', report)
    print(dict(cases=len(records), source=data['provenance']['source_sha256'], time_comparisons=times, controls=controls), flush=True)


if __name__ == '__main__':
    run()
