"""Scientific figure from complete replayed workspace and field-refinement runs."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Fresh report path required')
    inputs = {}
    def read(name):
        p = Path(name); value = json.loads(p.read_text()); inputs[str(p)] = sha(p)
        if not value.get('terminal'): raise ValueError('Complete input required: '+str(p))
        return value
    labels = ('Tensor', 'Release\nworkspace', 'Detach\ncache model', 'Axial\n288 cells')
    controls = (
        ('outputs/field_component_diagnostics/tensor_baseline_control/summary.json',
         'outputs/field_component_diagnostics/tensor_baseline_replay_control/summary.json'),
        ('outputs/workspace_memory_diagnostics/fine_control/summary.json',
         'outputs/workspace_memory_diagnostics/fine_replay_control/summary.json'),
        ('outputs/workspace_memory_diagnostics/detached_fine_control/summary.json',
         'outputs/workspace_memory_diagnostics/detached_fine_replay_control/summary.json'),
        ('outputs/workspace_memory_diagnostics/axial288_control/summary.json',
         'outputs/workspace_memory_diagnostics/axial288_replay_control/summary.json'))
    memory = []
    for pair in controls:
        rows = [read(name) for name in pair]
        if any(r['resource_guard_triggered'] or not r['all_cases_complete'] for r in rows):
            raise ValueError('Complete untriggered physical workflows required')
        if any(r['peak_rss_bytes'] > r['maximum_rss_bytes'] or
               r['minimum_free_bytes'] < r['threshold_free_bytes'] for r in rows):
            raise ValueError('Recorded resource limits must also hold at final observation')
        if not rows[1].get('complete_replay_verified'): raise ValueError('Full simulation replay required')
        verification = read(rows[1]['destination']+'_verification/summary.json')
        if not verification['bitwise_metrics']: raise ValueError('Bitwise full replay required')
        memory.append([r['peak_rss_bytes']/2**30 for r in rows])
    errors = []; strong = []
    for name in ('detached_baseline_curl', 'detached_axial288_curl'):
        s = read('outputs/'+name+'/summary.json')
        v = read('outputs/'+name+'_verification/summary.json')
        if not v['bitwise_metrics']: raise ValueError('Complete physical audit replay required')
        q = s['cases'][0]['quadrature'][-1]
        if q['order'] != 20 or not q['original_fingerprints_reproduced']:
            raise ValueError('Independent complete original-field reconstruction required')
        errors.append([100*q['metrics'][a]['relative_L2']['current_error'] for a in ('radial', 'azimuthal', 'axial')])
        strong.append([100*q['metrics'][a]['relative_L2']['strong_error'] for a in ('radial', 'azimuthal', 'axial')])
    identity = provenance(['examples.workspace_refinement_report'])
    destination.mkdir(parents=True); freeze(destination, identity)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))
    fig.subplots_adjust(left=.075, right=.97, bottom=.23, top=.77, wspace=.30)
    fig.suptitle('Workspace verification and axial field refinement', x=.075, ha='left', y=.965, fontsize=20)
    fig.text(.075, .89, 'Fixed fluid 64 / 256  •  radial field 256  •  independent radial-16, N128 reference', fontsize=11)
    colors = ('#225a80', '#168272'); width = .34; x = np.arange(4)
    m = np.asarray(memory)
    for k, label in enumerate(('Primary run', 'Full replay')):
        bars = axes[0].bar(x+(k-.5)*width, m[:, k], width, color=colors[k], label=label)
        axes[0].bar_label(bars, fmt='%.2f', padding=3, fontsize=9)
    axes[0].axhline(8, color='#a64635', ls='--', lw=1.2, label='Worker guard')
    axes[0].set(xticks=x, xticklabels=labels, ylim=(0, 9.5), ylabel='Own-process peak RSS (GiB)', title='Complete workflow memory')
    axes[0].legend(frameon=False, fontsize=9, loc='upper left')
    axes[0].grid(axis='y', alpha=.18); axes[0].set_axisbelow(True)
    x = np.arange(3)
    for k, label in enumerate(('Axial 256', 'Axial 288')):
        positions = x+(k-.5)*width
        axes[1].bar(positions, errors[k], width, color=colors[k], label=label)
        axes[1].scatter(positions, strong[k], s=35, facecolors='white', edgecolors='black', zorder=3,
                        label='Strong curl of magnetic rate' if k == 0 else None)
    axes[1].set(xticks=x, xticklabels=(r'$J_r$', r'$J_\theta$', r'$J_z$'), yscale='log',
                ylabel='Relative physical L2 error (%)', title='Each current-rate component')
    axes[1].legend(frameon=False, fontsize=9, loc='upper right')
    axes[1].grid(axis='y', which='both', alpha=.18); axes[1].set_axisbelow(True)
    fig.text(.075, .105, 'All full runs and physical audits replay bitwise; original numerical acceptance bounds remain.', fontsize=10)
    fig.text(.075, .065, 'Peak memory varies between runs. Strong curl excludes face distributions. Bounded periodic fixtures are not formation experiments.', fontsize=9, color='#444444')
    fig.savefig(destination/'verification.png', dpi=160)
    fig.savefig(destination/'verification.pdf'); plt.close(fig)
    if provenance(['examples.workspace_refinement_report']) != identity: raise ValueError('Source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], input_hashes=inputs,
        memory_GiB=memory, current_error_percent=errors, strong_curl_error_percent=strong,
        artifact_hashes={name: sha(destination/name) for name in ('verification.png', 'verification.pdf')},
        scope='Recorded full-workflow resource observations and independent q20 physical errors. Two field resolutions do not establish convergence order, time accuracy or formation.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
