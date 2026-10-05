"""Measured implicit wave convergence and unresolved contact-growth report."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .regular_plasma_reference import initialize
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    identity = provenance(['examples.implicit_plasma_report'])
    destination.mkdir(parents=True); freeze(destination, identity)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    colors = dict(midpoint='#3274a1', gauss='#33815b', radau='#c07527')
    summary = {}; inputs = {}
    for method, counts, template in [('midpoint', (30, 60, 120), 'implicit_wave_midpoint_{}'),
        ('gauss', (30, 60, 120, 240), 'implicit_wave_gauss_{}'), ('radau', (30, 60, 120, 240), 'radau_wave_{}')]:
        fields = []; records = []; audits = []
        for steps in counts:
            directory = Path('outputs')/template.format(steps)
            record = json.loads((directory/'case.json').read_text())
            audit_path = directory.with_name(directory.name+'_audit')/'summary.json'
            audit = json.loads(audit_path.read_text())
            raw = directory/record['arrays_file']
            if sha(raw) != record['arrays_sha256']: raise ValueError('Original wave arrays required')
            if not audit['full_replay'] or not audit['bitwise_metrics'] or not all(audit['bitwise_arrays'].values()):
                raise ValueError('Complete bitwise trajectory replay required')
            inputs[str(directory/'case.json')] = sha(directory/'case.json'); inputs[str(audit_path)] = sha(audit_path)
            arrays = dict(np.load(raw)); ref, d, _, system, _ = initialize(record['spec'])
            sample = TransportSampler(d, system.geometry.mesh, 16)
            field = sample.magnetic(arrays['final_magnetic']); fields.append(field)
            if len(fields) == 1:
                initial = sample.magnetic(arrays['initial_magnetic'])
                denominator = np.sqrt(np.sum(sample.weight*np.sum(initial*initial, axis=-1)))
            records.append(record); audits.append(audit)
        differences = [float(np.sqrt(np.sum(sample.weight*np.sum((a-b)**2, axis=-1)))/denominator)
                       for a, b in zip(fields, fields[1:])]
        ratios = (np.array(differences[:-1])/differences[1:]).tolist()
        summary[method] = dict(steps=counts, physical_B_time_differences=differences, time_ratios=ratios,
            physical_energy_over_wave_energy=[a['order20']['change_over_wave_energy'] for a in audits],
            grid_energy_over_wave_energy=[a['energy_change_over_wave_energy'] for a in records],
            final_B_relative_errors=[a['final']['magnetic']['relative_L2'] for a in records])
        label=f'{method.capitalize()} (ratio {ratios[-1]:.3f})'
        axes[0, 0].loglog(counts[:-1], differences, 'o-', color=colors[method], label=label)
        axes[0, 1].loglog(counts, np.abs(summary[method]['physical_energy_over_wave_energy']), 'o-', color=colors[method], label=method.capitalize())
        axes[0, 2].loglog(counts, summary[method]['final_B_relative_errors'], 'o-', color=colors[method], label=method.capitalize())
    axes[0, 0].set(title='Time refinement: full wave periods', xlabel='Steps in coarser run', ylabel='Successive physical B difference / initial B norm')
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set(title='Physical energy drift, separate quadrature', xlabel='Steps per period', ylabel='|Energy change| / wave energy')
    axes[0, 1].legend(fontsize=8)
    axes[0, 2].set(title='Total B error against continuum wave', xlabel='Steps per period', ylabel='Relative physical L2 error')
    axes[0, 2].legend(fontsize=8)
    contacts = [('implicit_gauss_newton_contact', 'Gauss: 4 elements, 120 steps'),
                ('implicit_gauss_contact_time240', 'Gauss: 4 elements, 240 steps'),
                ('implicit_gauss_contact_ne8', 'Gauss: 8 elements, 120 steps'),
                ('radau_contact_ne4_t120', 'Radau: 4 elements, 120 steps')]
    # Only the last saved checkpoint is retained by these stopped pilots.
    for name, label in contacts:
        path = Path('outputs')/name/'termination.json'; record = json.loads(path.read_text())
        inputs[str(path)] = sha(path); last = record['saved_checkpoint']; t = last['time']*1e6
        axes[1, 0].plot(t, last['density_min']/1e20, 'o', label=label)
        axes[1, 1].semilogy(t, last['grid_change']['magnetic'], 'o', label=label)
    axes[1, 0].axhline(.75, color='black', linestyle='--', label='Exact continuous minimum')
    axes[1, 0].set(title='Stopped contacts: last saved density', xlabel='Saved time (microseconds)', ylabel='Minimum point density / reference density', ylim=(.3, .8))
    axes[1, 0].legend(fontsize=7, loc='lower left')
    axes[1, 1].set(title='Stopped contacts: unwanted magnetic energy', xlabel='Saved time (microseconds)', ylabel='Incremental grid magnetic energy (J)')
    axes[1, 1].text(.04, .06, 'Exact contact: zero increment\nNo run completed the 100 microsecond transit', transform=axes[1, 1].transAxes, fontsize=8)
    path = Path('outputs/stationary_contact_audit/summary.json'); spectrum = json.loads(path.read_text())
    inputs[str(path)] = sha(path); subset = [v for v in spectrum['continuum'] if v['radial_mode'] == 1]
    axes[1, 2].semilogx([v['points'] for v in subset], [v['maximum_real']/1000 for v in subset], 'o-', color='#715595')
    axes[1, 2].set(title='Separate continuum contact linearization', xlabel='Axial points (first Bessel radial mode)', ylabel='Largest real growth rate (thousand / s)')
    axes[1, 2].text(.04, .92, 'Different radial spaces from production probe;\nnot a matched-mode comparison', va='top', transform=axes[1, 2].transAxes, fontsize=8)
    for ax in axes.ravel(): ax.grid(alpha=.2)
    fig.suptitle('Verified time integration; long-contact stability remains unresolved', fontsize=17, weight='bold')
    fig.savefig(destination/'convergence.png', dpi=150); fig.savefig(destination/'convergence.pdf'); plt.close(fig)
    if provenance(['examples.implicit_plasma_report']) != identity: raise RuntimeError('Report source changed')
    (destination/'summary.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'], inputs=inputs,
        terminal=True, waves=summary, tests=dict(passed=432, seconds=227.47),
        limitations='Filled periodic verification. No completed accurate contact transit, exact vacuum, formation, or laboratory design validation.'), indent=2)+'\n')


if __name__ == '__main__': run('outputs/implicit_plasma_report')
