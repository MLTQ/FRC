"""Scientific report for exact curl actions and complete physical mode spectra."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, curl, spectra):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    inputs = {}; cases = {}
    for directory in spectra:
        root = Path(directory); record = json.loads((root/'summary.json').read_text())
        if not record['terminal']: raise ValueError('Complete spectra required')
        if sha(root/record['arrays_file']) != record['arrays_sha256']: raise ValueError('Unchanged spectrum arrays required')
        inputs[str(root/'summary.json')] = sha(root/'summary.json')
        inputs[str(root/record['arrays_file'])] = record['arrays_sha256']
        data = dict(np.load(root/record['arrays_file'])); modes = data['signed_cell_mode']
        q = next(v for v in record['quadrature'] if v['order'] == 20)
        for amplitude, key in zip(record['amplitudes'], ('first_error', 'second_error')):
            item = dict(spectrum=data[f'q20_current_radial_{key}_spectrum'],
                reference=data['q20_current_radial_reference_spectrum'],
                integral=q['metrics']['current_radial'][key]['direct_squared_integral'],
                high_fraction=q['metrics']['current_radial'][key]['fractions']['above_reference_band'])
            if amplitude in cases and any(not np.array_equal(item[k], cases[amplitude][k]) for k in item):
                raise ValueError('Same amplitude must reproduce the same physical spectrum')
            cases[amplitude] = item
    root = Path(curl); record = json.loads((root/'summary.json').read_text())
    if not record['terminal']: raise ValueError('Complete curl action comparison required')
    inputs[str(root/'summary.json')] = sha(root/'summary.json')
    result = record['observations'][-1]; identity = provenance(['examples.curl_spectrum_report'])
    destination.mkdir(parents=True); freeze(destination, identity)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    methods = ['assembled', 'dense_factored', 'paired_factored', 'exact_raw_factored']
    labels = ['Assembled', 'Factored', 'Paired stencil', 'Exact raw curl']
    x = np.arange(len(methods)); a = axes[0, 0]
    for offset, key, label in [(-.18, 'distribution_over_rhs', 'Electric distributivity'),
        (.18, 'distribution_plus_low_over_rhs', 'After addition-low accounting')]:
        a.bar(x+offset, [result['metrics'][k][key] for k in methods], .36, label=label)
    a.set(yscale='log', xticks=x, xticklabels=labels, ylabel='Residual norm / original RHS norm',
        title='A  Same rejected iterate, alternate arithmetic', ylim=(1e-16, 3e-10))
    a.tick_params(axis='x', labelrotation=15); a.legend(fontsize=8, loc='lower left')
    a = axes[0, 1]
    for amplitude, values in sorted(cases.items(), reverse=True):
        folded = np.array([np.sum(values['spectrum'][abs(modes) == k]) for k in range(max(abs(modes))+1)])
        a.semilogy(np.arange(len(folded)), np.sqrt(folded/np.sum(values['reference'])), label=f'epsilon = {amplitude:g}')
    a.axvline(64, color='0.5', ls=':', label='Reference band edge')
    a.set(xlabel='Absolute axial cell-mode index', ylabel='Radial-current error norm / reference norm',
        title='B  Full physical radial-current error spectrum', ylim=(1e-9, .1))
    a.legend(fontsize=8)
    a = axes[1, 0]; amplitudes = sorted(cases, reverse=True); positions = np.arange(len(amplitudes))
    total = np.array([np.sqrt(cases[v]['integral']/np.sum(cases[v]['reference'])) for v in amplitudes])
    fraction = np.array([cases[v]['high_fraction'] for v in amplitudes])
    a.bar(positions-.22, total*100, .22, label='Total')
    a.bar(positions, total*np.sqrt(1-fraction)*100, .22, label='Within reference band')
    a.bar(positions+.22, total*np.sqrt(fraction)*100, .22, label='Above reference band')
    a.set(xticks=positions, xticklabels=[f'{v:g}' for v in amplitudes], xlabel='Centered perturbation amplitude',
        ylabel='Relative radial-current error (%)', title='C  Component accuracy remains unresolved')
    a.legend(fontsize=8)
    a = axes[1, 1]; errors = result['raw_error_against_exact_over_rhs']; x = np.arange(2)
    for offset, key, label in [(-.18, 'raw', 'Original CSR raw curl'), (.18, 'paired_raw', 'Subtraction-first raw curl')]:
        a.bar(x+offset, [errors[v][key] for v in ('e0', 'e')], .36, label=label)
    a.set(yscale='log', xticks=x, xticklabels=['Forcing electric term', 'Rounded full electric field'],
        ylabel='Projected raw-curl error / original RHS norm', title='D  Independent exact raw-product comparison')
    a.legend(fontsize=8, loc='lower left')
    for a in axes.flat: a.grid(axis='y', alpha=.2); a.set_axisbelow(True)
    fig.suptitle('FRC numerical verification: curl arithmetic and individual current accuracy', fontsize=14)
    fig.savefig(destination/'verification.png', dpi=160); plt.close(fig)
    if provenance(['examples.curl_spectrum_report']) != identity: raise ValueError('Source changed')
    (destination/'summary.json').write_text(json.dumps(dict(terminal=True, source_sha256=identity['source_sha256'],
        input_hashes=inputs, amplitudes=amplitudes, figure_sha256=sha(destination/'verification.png'),
        scope='Filled-periodic verification fixture; no solver repair, formation or confinement claim'), indent=2)+'\n')
