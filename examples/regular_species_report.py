"""Measured species convergence figure and machine-readable checkpoint."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def report(source, output):
    source, output = Path(source), Path(output)
    audit = json.loads((source/'audit.json').read_text())
    identity = json.loads((source/'provenance.json').read_text())
    if not audit['terminal'] or audit['source_sha256'] != identity['source_sha256']:
        raise ValueError('Complete matching audit required')
    records = [json.loads((source/row['file']).read_text()) for row in audit['records']]
    for record in records:
        raw = source/record['arrays_file']
        assert hashlib.sha256(raw.read_bytes()).hexdigest() == record['arrays_sha256']
    output.mkdir(parents=True, exist_ok=True)
    if (output/'convergence.png').exists() or (output/'summary.json').exists():
        raise ValueError('Report already exists')
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.7), constrained_layout=True)
    colors = ('#126D9A', '#BD5C17')
    for p, style in ((3, '--'), (5, '-')):
        ids = [i for i, r in enumerate(records) if r['kind'] == 'rate' and r['spec'] == dict(elements=r['spec']['elements'], degree=p)]
        ns = [records[i]['spec']['elements'] for i in ids]
        for k, name in enumerate(('electron', 'ion')):
            for ax, metric in ((axes[0, 0], 'entropy'), (axes[0, 1], 'pressure')):
                ys = [audit['records'][i]['order20']['species'][k][metric]['relative_L2'] for i in ids]
                ax.loglog(ns, ys, style+'o', color=colors[k], label=f'{name}, degree {p}')
    ids = list(range(4, 8))
    ns = [records[i]['spec']['elements'] for i in ids]
    for name, color in zip(('axis_element', 'wall_element', 'maximum_over_rms'), ('#126D9A', '#BD5C17', '#387E51')):
        axes[0, 2].loglog(ns, [audit['records'][i]['order20']['species'][0]['entropy'][name] for i in ids], '-o', color=color, label=name.replace('_', ' '))
    ids = (14, 15, 16)
    ns = [records[i]['spec']['elements'] for i in ids]
    for name, color in zip(('relative_L2', 'axis_element', 'wall_element'), ('#126D9A', '#BD5C17', '#387E51')):
        axes[1, 0].loglog(ns, [audit['records'][i]['order20']['signal_error'][name] for i in ids], '-o', color=color, label=name.replace('_', ' '))
    ids = list(range(4, 8))
    ns = [records[i]['spec']['elements'] for i in ids]
    axes[1, 1].loglog(ns, [records[i]['controls']['point_skew'] for i in ids], '--o', color='#B84045', label='point transfer + skew')
    axes[1, 1].loglog(ns, [records[i]['species'][0]['point_rate']['relative_L2'] for i in ids], '-o', color='#126D9A', label='physical projection + skew')
    ids = (14, 15, 16)
    ns = [records[i]['spec']['elements'] for i in ids]
    axes[1, 2].loglog(ns, [abs(records[i]['grid_quadratic_relative_drift']) for i in ids], '-o', color='#126D9A', label='grid quadratic drift')
    axes[1, 2].loglog(ns, [abs(audit['records'][i]['order20']['physical_quadratic_relative_drift']) for i in ids], '-o', color='#BD5C17', label='continuous quadratic drift')
    titles = ('Original entropy-rate error', 'Original pressure-field error', 'Electron entropy: local errors',
              'One full sheared transit: signal error', 'Retained inaccurate transfer control', 'Grid and continuous measures differ')
    for ax, title in zip(axes.ravel(), titles):
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Radial fluid elements')
        ax.set_ylabel('Relative error / absolute relative drift')
        ax.grid(True, which='both', alpha=.18)
        ax.legend(fontsize=8)
    fig.suptitle('Regular species coupling | bounded verification\nIndependent equations, fixed-current evolution and explicit accuracy limits', fontsize=15)
    fig.savefig(output/'convergence.png', dpi=180)
    fig.savefig(output/'convergence.pdf')
    plt.close(fig)
    steps = [17, 18, 15, 19]
    endpoints = [dict(np.load(source/records[i]['arrays_file']))['final_entropy'] for i in steps]
    differences = [float(np.linalg.norm(b-a)) for a, b in zip(endpoints[:-1], endpoints[1:])]
    time_ratios = [a/b for a, b in zip(differences[:-1], differences[1:])]
    changes = []
    for row in audit['records']:
        if row['kind'] != 'rate':
            continue
        for k in range(2):
            for metric in ('entropy', 'pressure'):
                low, high = row['order16']['species'][k][metric], row['order20']['species'][k][metric]
                # Maxima use different sample locations and are deliberately excluded.
                for name in ('relative_L2', 'axis_element', 'wall_element'):
                    changes.append(abs(high[name]-low[name])/max(abs(high[name]), abs(low[name])))
    from .regular_species_reference import initialize
    _, fluid, state, _, _ = initialize(records[23]['spec'])
    root = np.sqrt(fluid.volume*state[..., 0]).ravel()
    control_arrays = dict(np.load(source/records[23]['arrays_file']))
    constant_residuals = {}
    for label in records[23]['spectra']:
        matrix = control_arrays[label+'_operator']
        constant_residuals[label] = float(np.linalg.norm(matrix@root)/(np.linalg.norm(matrix)*np.linalg.norm(root)))
    summary = dict(source_sha256=identity['source_sha256'], source_files=len(identity['source_files']),
                   cases=len(records), rate_cases=14, full_transit_cases=8, spectral_controls=2,
                   all_arrays_and_metrics_bitwise=all(r['arrays_and_metrics_bitwise'] for r in audit['records']),
                   time_steps=[records[i]['steps'] for i in steps], time_differences=differences, time_ratios=time_ratios,
                   maximum_relative_quadrature_change_in_integrated_error=float(max(changes)),
                   finest_rate=audit['records'][7]['order20'], finest_transit=audit['records'][16]['order20'],
                   max_work_backward=max(s['work_backward'] for r in records if r['kind']=='rate' for s in r['species']),
                   max_grid_power_backward=max(r['total_grid_power_backward'] for r in records if r['kind']=='rate'),
                   max_original_ohm_backward=max(r['original_ohm_backward'] for r in records if r['kind']=='rate'),
                   control_spectra=records[23]['spectra'],
                   constant_state_relative_residual=constant_residuals,
                   report_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('source')
    p.add_argument('output')
    args = p.parse_args()
    report(args.source, args.output)
