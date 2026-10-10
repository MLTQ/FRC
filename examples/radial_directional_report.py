"""Measured directional refinement, cancellation and strict arithmetic checks."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from .implicit_plasma_records import provenance, freeze, sha


def run(destination='outputs/radial16_directional_report'):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    names = ['radial16_balanced_rates', 'radial16_directional_refinement',
             'radial16_fluid_refinement', 'radial16_term_audit',
             'radial16_drive_audit', 'radial16_inverse_actions',
             'radial16_split_inverse_probe', 'radial16_split_inverse_rates']
    paths = {name: Path('outputs')/name/'summary.json' for name in names}
    data = {name: json.loads(path.read_text()) for name, path in paths.items()}
    if not all(record['terminal'] for record in data.values()):
        raise ValueError('Terminal input studies required')
    identity = provenance(['examples.radial_directional_report'])
    destination.mkdir(parents=True); freeze(destination, identity)
    cases = [data[names[0]]['cases'][1], data[names[1]]['cases'][0],
             data[names[2]]['cases'][0], data[names[7]]['cases'][1]]
    if any(case['status'] != 'complete' for case in cases):
        raise ValueError('Only accepted full rates may appear in rate comparison')
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 12,
                        'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(3, 2, figsize=(14, 13), layout='constrained')
    colors = ['#176B87', '#C85C32', '#617A35', '#825D99']
    labels = ['Baseline\nF32 B128 z128', 'Axial doubled\nF32 B128 z256',
              'Fluid doubled\nF64 B128 z128', 'Field radial doubled*\nF32 B256 z128']
    ax = axes[0, 0]; x = np.arange(4)
    for j, (key, label) in enumerate(zip(
            ['radial_velocity', 'radial_B', 'current_rate', 'azimuthal_B'],
            ['Radial velocity', 'Radial B', 'Current rate', 'Toroidal B'])):
        ax.bar(x+(j-1.5)*.18, [c['quadrature'][-1]['rate_errors'][key]['relative_L2']
                             for c in cases], width=.18, label=label, color=colors[j])
    ax.set_yscale('log'); ax.set_xticks(x, labels, fontsize=8)
    ax.set_ylabel('Relative physical L2 rate error'); ax.legend(fontsize=8)
    ax.set_title('A. Separate mesh controls expose remaining accuracy error')
    metric = data[names[3]]['cases'][1]['quadrature'][-1]['metrics']
    terms = ['ideal', 'hall', 'pressure_full', 'inertia']; labels2 = ['Ideal', 'Hall', 'Pressure', 'Inertia']
    net = metric['faraday_accepted_rate']['volume']['target_L2'][1]
    ax = axes[0, 1]
    ax.bar(range(4), [metric['faraday_'+t]['volume']['target_L2'][1]/net for t in terms], color=colors)
    ax.set_yscale('log'); ax.set_xticks(range(4), labels2)
    ax.set_ylabel('Term L2 / net toroidal-B rate L2')
    ax.set_title('B. Physical terms are much larger than their net response')
    ax = axes[1, 0]
    ax.bar(range(5), [metric['faraday_'+t]['volume']['error_over_net_L2'][1]
                     for t in terms+['accepted_rate']], color=colors+['#333333'])
    ax.set_yscale('log'); ax.set_xticks(range(5), labels2+['Accepted net'])
    ax.set_ylabel('Term error / net toroidal-B rate L2')
    ax.set_title('C. Term errors also cancel; inertia includes J-rate feedback')
    rawpath = Path('outputs/radial16_term_audit/case_01.npz')
    raw = dict(np.load(rawpath)); radius = raw['q20_radius']
    ax = axes[1, 1]
    for term, label, color in [('lorentz', 'Lorentz', colors[1]), ('bulk', 'Pressure', colors[0]),
                               ('accepted_rate', 'Complete momentum', '#333333')]:
        values = raw['q20_momentum_'+term][:, 2, 0]
        ax.plot(radius/(radius[-1]+radius[0]), values/np.sum(values),
                label=label, color=color)
    ax.axvspan(31/32, 1, color='#777777', alpha=.13, label='Outer fluid element')
    ax.set_xlabel('r / tube radius'); ax.set_ylabel('Fraction of squared radial-rate error per field cell')
    ax.set_title('D. Most radial momentum error is near the outer wall'); ax.legend(fontsize=8)
    ax = axes[2, 0]
    drive = data[names[4]]['cases']
    for term, label, color in [('hall', 'Hall', colors[1]), ('pressure_reduced', 'Pressure', colors[0]),
                               ('driving_sum', 'Sum of driving terms', '#333333')]:
        ax.semilogy(range(3), [c['quadrature'][-1]['metrics']['driving_'+term]['volume']['error_over_full_reference_rate_L2'][1]
                              for c in drive], 'o-', label=label, color=color)
    ax.set_xticks(range(3), ['F32 B128', 'F32 B256\noriginal solve rejected', 'F64 B128'])
    ax.set_ylabel('Driving error / full reference B-rate L2')
    ax.set_title('E. Driving-only audit never accepts a magnetic rate'); ax.legend(fontsize=8)
    ax = axes[2, 1]; records = data[names[5]]['records']
    for key, label, color in [('observed_gap_over_rhs', 'Observed residual gap', '#333333'),
                               ('mass_gap_over_rhs', 'Mass-inverse distributivity', colors[1]),
                               ('current_gap_over_rhs', 'Current association', colors[0])]:
        ax.semilogy(range(len(records)), [r[key] for r in records], 'o-', label=label, color=color)
    ax.axhline(2e-10, ls='--', color='#777777', label='Twice original dual bound')
    ax.set_xlabel('Saved original refinement iterate'); ax.set_ylabel('Gap / original RHS norm')
    ax.set_title('F. Independent reconstruction of the failed residual gap'); ax.legend(fontsize=8)
    for ax in axes.ravel(): ax.grid(axis='y', alpha=.2)
    fig.suptitle('Fixed radial16 verification: separate spatial error from roundoff', fontsize=17)
    fig.supxlabel('Filled periodic physical-mass fixture; F = fluid elements, B = radial field cells.\n'
                  '*Opt-in split-inverse arithmetic; all original bounds retained. Formation and full trajectory accuracy remain open.', fontsize=11)
    for extension in ('png', 'pdf'): fig.savefig(destination/f'verification.{extension}', dpi=150)
    plt.close(fig)
    inputs = {str(path): sha(path) for path in paths.values()}; inputs[str(rawpath)] = sha(rawpath)
    record = dict(terminal=True, source_sha256=identity['source_sha256'], inputs=inputs,
                  selected_specs=[case['spec'] for case in cases],
                  figures={f'verification.{e}': sha(destination/f'verification.{e}') for e in ('png', 'pdf')})
    (destination/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
