"""Scientific report of complete spatial workflows and trajectory integration."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, spatial, trajectory, energy):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing overwrite')
    hashes = {}
    def load(directory):
        path = Path(directory)/'summary.json'
        value = json.loads(path.read_text())
        if not value['terminal']:
            raise ValueError('Terminal scientific input required')
        hashes[str(path)] = sha(path)
        return value
    rows = []
    for directory, index, label in spatial:
        case = load(directory)['cases'][index]
        if case['status'] != 'complete':
            raise ValueError('Complete signed physical case required')
        raw = Path(directory)/case['arrays_file']
        if sha(raw) != case['arrays_sha256']:
            raise ValueError('Original full physical arrays required')
        hashes[str(raw)] = sha(raw)
        rows.append((label, case['quadrature'][-1]))
    audit = load(trajectory); scale = load(energy)['cases'][0]
    if scale['input_summary_sha256'] != sha(Path(trajectory)/'summary.json'):
        raise ValueError('Energy normalization must belong to this trajectory audit')
    identity = provenance(['examples.joint_e2e_report'])
    destination.mkdir(parents=True); freeze(destination, identity)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.titleweight': 'bold', 'axes.grid': True, 'grid.alpha': .18})
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    colors = ['#2864b4', '#db6e28', '#2b9178', '#8c56a3']
    x = np.arange(len(rows))
    for name, label, color in zip(('radial_velocity', 'radial_B', 'azimuthal_B', 'current_rate'),
                                  ('Radial velocity', 'Radial B', 'Toroidal B', 'Current'), colors):
        axes[0, 0].semilogy(x, [q['rate_errors'][name]['relative_L2'] for _, q in rows], 'o-', label=label, color=color)
    axes[0, 0].set(title='Fixed radial16: complete rate error', ylabel='Relative physical L2 error')
    axes[0, 0].set_xticks(x, [label for label, _ in rows]); axes[0, 0].legend(fontsize=9)
    for name, label, color in [('initial_current', 'Initial weak current', colors[0]),
                                ('strong_ampere', 'Weak versus strong curl', colors[1])]:
        axes[0, 1].semilogy(x, [q['representation_and_faces']['volume'][name]['relative_L2'] for _, q in rows],
                           'o-', label=label, color=color)
    axes[0, 1].semilogy(x, [q['representation_and_faces']['internal_faces']['radial']['normal_current']['jump_over_reference_L2']
                           for _, q in rows], 'o-', label='Internal normal-current jump', color=colors[2])
    axes[0, 1].set(title='Initial current representation and traces', ylabel='Relative L2 error / jump')
    axes[0, 1].set_xticks(x, [label for label, _ in rows]); axes[0, 1].legend(fontsize=9)
    final = audit['quadrature'][-1]['physical']['final']
    names = ['density', 'radial_velocity', 'radial_B', 'azimuthal_B', 'current', 'electric']
    values = [final['errors'][name]['relative_L2'] for name in names]
    axes[1, 0].bar(np.arange(len(names)), values, color=colors[0], width=.6)
    axes[1, 0].set_yscale('log')
    axes[1, 0].set_xticks(np.arange(len(names)), ['Density', 'Radial v', 'Radial B', 'Toroidal B', 'Current', 'Electric'])
    axes[1, 0].set(title='Separate radial1 trajectory: final error', ylabel='Relative physical L2 error', ylim=(1e-4, 3e-2))
    for k, value in enumerate(values):
        axes[1, 0].text(k, value*1.12, f'{100*value:.3f}%', ha='center', fontsize=9)
    for q, marker in zip(scale['quadrature'], ('o', 'x')):
        values = [abs(q['signed_energy_drift_over_initial_quadratic'][name]) for name in ('minus', 'plus')]
        axes[1, 1].semilogy([0, 1], values, marker, ms=8, label=f"Quadrature {q['order']}")
    axes[1, 1].set_xticks([0, 1], ['Negative perturbation', 'Positive perturbation'])
    axes[1, 1].set(title='Actual total-energy drift / initial signal scale', ylabel='Absolute drift ratio')
    axes[1, 1].set_xlim(-.4, 1.4); axes[1, 1].legend(fontsize=9)
    fig.suptitle('FRC numerical verification — combined spatial refinement and full evolution', fontsize=15, fontweight='bold')
    fig.text(.5, .934, 'Top: fixed radial16,N128 reference. F = fluid elements/axial points; B = field cells/axial points.', ha='center', fontsize=10)
    fig.text(.5, .032, 'Bottom: both signed 100 μs trajectories, 60 Gauss steps, F4/32 and B8/32. Signal energy is a scale, not a conserved invariant.', ha='center', fontsize=9)
    fig.text(.5, .012, 'Filled periodic verification fixtures. Mesh/time convergence, vacuum coupling and FRC formation remain open.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .055, 1, .917), h_pad=2.5, w_pad=2.5)
    for suffix in ('png', 'pdf'):
        fig.savefig(destination/f'verification.{suffix}', dpi=180)
    plt.close(fig)
    result = dict(terminal=True, source_sha256=identity['source_sha256'], input_hashes=hashes,
        spatial_inputs=spatial, trajectory_audit=trajectory, energy_scale=energy,
        scope='Instantaneous radial16 spatial consistency and separate full radial1 integration; not formation')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
