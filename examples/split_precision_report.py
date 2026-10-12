"""Scientific view of physical component errors and residual arithmetic."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .finite_mass_contact_reference import FiniteMassContact
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, components='outputs/component_amplitude_localization',
        arithmetic='outputs/split_residual_arithmetic_audit', probe='outputs/split_residual_probe'):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    inputs = {}; records = []
    for directory in (components, arithmetic, probe):
        root = Path(directory); path = root/'summary.json'
        record = json.loads(path.read_text()); source = json.loads((root/'provenance.json').read_text())
        if not record['terminal'] or record['source_sha256'] != source['source_sha256']:
            raise ValueError('Terminal source-identified input required')
        if any(sha(p) != h for p, h in source['source_files'].items()): raise ValueError('Original source changed')
        if sha(root/record['arrays_file']) != record['arrays_sha256']: raise ValueError('Input array hash mismatch')
        inputs[str(path)] = sha(path); inputs[str(root/record['arrays_file'])] = record['arrays_sha256']
        records.append(record)
    component, exact, observed = records
    if not exact['exact_electric_sum_verified'] or not observed['observation_consistent']:
        raise ValueError('Matched original observation and independent arithmetic audit required')
    q = component['quadrature'][-1]; metrics = q['metrics']
    data = dict(np.load(Path(components)/component['arrays_file']))
    identity = provenance(['examples.split_precision_report']); destination.mkdir(parents=True); freeze(destination, identity)
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 10})
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    colors = ('#2166ac', '#d6604d'); labels = ('amplitude 1e-4', 'amplitude 1e-5')
    names = ('current_radial', 'current_azimuthal', 'current_axial', 'current_rate')
    y = np.arange(len(names)); ax = axes[0, 0]
    for index, key in enumerate(('first_error', 'second_error')):
        values = [100*metrics[n]['regions']['volume']['relative_L2'][key] for n in names]
        ax.barh(y+(index-.5)*.3, values, height=.3, color=colors[index], label=labels[index])
    ax.set_yticks(y, ('radial J rate', 'azimuthal J rate', 'axial J rate', 'full J-rate vector'))
    ax.invert_yaxis(); ax.set_xscale('log'); ax.set_xlim(.01, 12)
    ax.set_xlabel('Physical L2 error / own continuum reference (%)')
    ax.set_title('A. The total norm hides component error'); ax.legend(loc='lower right', fontsize=9)
    ax.grid(axis='x', alpha=.25)
    ax = axes[0, 1]; names2 = ('current_radial', 'current_axial', 'radial_velocity', 'azimuthal_B')
    bottom = np.zeros(4)
    for region, color, label in [('axis', '#8dd3c7', 'axis: inner 5% R'), ('interior', '#80b1d3', 'interior'), ('wall', '#fb8072', 'wall: outer 5% R')]:
        values = np.array([100*metrics[n]['regions'][region]['fraction_of_domain_squared']['difference'] for n in names2])
        ax.bar(np.arange(4), values, bottom=bottom, color=color, label=label); bottom += values
    ax.set_xticks(range(4), ('radial J rate', 'axial J rate', 'radial velocity', 'toroidal B rate'))
    ax.set_ylim(0, 110); ax.set_ylabel('Share of squared amplitude difference (%)')
    ax.set_title('B. Sensitivity mainly occupies the interior'); ax.legend(loc='lower left', fontsize=9)
    ax = axes[1, 0]; order = q['order']; radius = FiniteMassContact().radius
    r = data[f'q{order}_radius'].reshape(-1, order).mean(axis=1)/radius
    for key, color, label in [('first_error', colors[0], 'error, amplitude 1e-4'), ('second_error', colors[1], 'error, amplitude 1e-5'), ('difference', '#444444', 'amplitude difference')]:
        values = data[f'q{order}_current_radial_{key}_radial_integrals'].reshape(-1, order).sum(axis=1)
        ax.plot(r, values/values.sum()*len(r), color=color, label=label, lw=1.1)
    ax.axvspan(0, .05, color='#8dd3c7', alpha=.25); ax.axvspan(.95, 1, color='#fb8072', alpha=.2)
    ax.set_xlim(0, 1); ax.set_xlabel('Radius / chamber radius')
    ax.set_ylabel('Normalized squared-error density per r/R')
    ax.set_title('C. Radial-current errors extend across the radius'); ax.legend(fontsize=8)
    ax = axes[1, 1]; last = exact['observations'][-1]; obs = observed['observations'][-1]
    values = [last['gap_over_rhs'], last['terms']['electric_curl_distribution']['norm_over_rhs'],
        obs['electric_addition_prediction_over_rhs'], last['terms']['inertia_paths']['norm_over_rhs'],
        last['terms']['metric_association']['norm_over_rhs'], obs['decomposition_remainder_over_rhs']]
    text = ('observed vector gap', 'electric curl distribution', 'electric addition only',
        'current/inertia paths', 'H/T metric association', 'all-term remainder')
    ax.barh(np.arange(6), values, color=('#333333', '#b2182b', '#ef8a62', '#67a9cf', '#2166ac', '#999999'))
    ax.set_yticks(np.arange(6), text); ax.invert_yaxis(); ax.set_xscale('log')
    ax.set_xlim(1e-17, 2e-10); ax.set_xlabel('Vector norm / original right-hand-side norm')
    ax.set_title('D. Final rejected iterate: separate diagnostics'); ax.grid(axis='x', alpha=.25)
    fig.suptitle('Fixed radial16 mode: component accuracy and residual arithmetic', fontsize=16, y=.98)
    fig.text(.5, .942, 'Filled periodic verification fixture; no formation or confinement claim', ha='center', color='#444444')
    fig.text(.5, .018, 'Electric-addition loss alone leaves 66.2% of the gap. Regional fractions use masked original quadrature nodes.', ha='center', fontsize=10)
    fig.subplots_adjust(left=.12, right=.98, top=.89, bottom=.09, wspace=.6, hspace=.48)
    fig.savefig(destination/'verification.png', dpi=180); fig.savefig(destination/'verification.pdf'); plt.close(fig)
    result = dict(terminal=True, source_sha256=identity['source_sha256'], input_hashes=inputs,
        scope='Scientific rendering of completed observations; simulation replay records are separate',
        png_sha256=sha(destination/'verification.png'), pdf_sha256=sha(destination/'verification.pdf'))
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
