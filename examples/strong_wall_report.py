"""Measured wall correction, full contacts and remaining continuum-rate limits."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .strong_wall_reference import initialize
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance, freeze, sha


def load_case(kind, elements, steps, inputs):
    root = Path(f'outputs/conforming_{kind}_ne{elements}_t{steps}')
    case = json.loads((root/'case.json').read_text())
    path = root.with_name(root.name+'_audit')/'summary.json'; audit = json.loads(path.read_text())
    raw = root/case['arrays_file']
    if sha(raw) != case['arrays_sha256']: raise ValueError('Original raw hash required')
    if not audit['terminal'] or not audit['full_replay'] or not audit['bitwise_metrics'] or not all(audit['bitwise_arrays'].values()):
        raise ValueError('Complete physical audit and bitwise replay required')
    inputs[str(path)] = sha(path); inputs[str(root/'case.json')] = sha(root/'case.json')
    return case, audit, dict(np.load(raw))


def temporal(kind, elements, counts, inputs):
    values = [load_case(kind, elements, steps, inputs) for steps in counts]
    ref, d, _, system, _ = initialize(values[0][0]['spec']); records = {}
    for order in (16, 20):
        sample = TransportSampler(d, system.geometry.mesh, order)
        fields = [sample.magnetic(a['final_magnetic']) if kind == 'wave' else
                  sample.scalar(a['final_state'][..., 0]-a['initial_state'][..., 0]) for _, _, a in values]
        if kind == 'wave':
            initial = sample.magnetic(values[0][2]['initial_magnetic'])
            denominator = np.sqrt(np.sum(sample.weight*np.sum(initial**2, axis=-1)))
            square = lambda x: np.sum(x*x, axis=-1)
        else:
            denominator = np.sqrt(np.sum(sample.weight*(ref.fields(sample.Z)['number_density']-ref.number_density)**2))
            square = lambda x: x*x
        differences = [float(np.sqrt(np.sum(sample.weight*square(a-b)))/denominator) for a, b in zip(fields, fields[1:])]
        records[str(order)] = dict(differences=differences, ratios=(np.array(differences[:-1])/differences[1:]).tolist())
    return dict(steps=counts, physical_time_differences=records,
        continuous_energy_relative=[a['order20']['relative_energy_drift'] for _, a, _ in values],
        grid_energy_relative=[r['grid_energy_relative_drift'] for r, _, _ in values])


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    identity = provenance(['examples.strong_wall_report']); inputs = {}
    wall_path = Path('outputs/contact_wall_diagnosis_conforming/summary.json')
    wall = json.loads(wall_path.read_text()); inputs[str(wall_path)] = sha(wall_path)
    mode_path = Path('outputs/continuum_wall_probe/summary.json')
    mode = json.loads(mode_path.read_text()); inputs[str(mode_path)] = sha(mode_path)
    for name in ('contact_wall_replay_verification', 'continuum_wall_replay_verification'):
        path = Path('outputs')/name/'summary.json'
        repeated = json.loads(path.read_text()); inputs[str(path)] = sha(path)
        if not repeated['terminal'] or not repeated['bitwise_metrics'] or not all(v['bitwise_equal'] for v in repeated['arrays'].values()):
            raise ValueError('Complete diagnostic replay required')
    contacts = [load_case('contact', n, 120, inputs) for n in (2, 4, 8)]
    waves = [load_case('wave', n, 120, inputs) for n in (2, 4, 8)]
    time = dict(wave=temporal('wave', 2, [30, 60, 120, 240], inputs),
                contact=temporal('contact', 4, [60, 120, 240], inputs))
    destination.mkdir(parents=True); freeze(destination, identity)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    plt.rcParams.update({'font.size': 10})
    raw = wall_path.parent/wall['physical_mode'][0]['arrays_file']
    if sha(raw) != wall['physical_mode'][0]['arrays_sha256']: raise ValueError('Original mode reconstruction required')
    a = dict(np.load(raw)); radial = a['radius']; physical = a['physical_mode']
    for k, label in [(0, 'Density'), (1, 'Radial velocity'), (7, 'Toroidal B')]:
        profile = np.sum(a['axial_weight'][None, :]*physical[..., k]**2, axis=1)
        axes[0, 0].plot(radial/.035, profile/profile.max(), label=label)
    axes[0, 0].set(title='Original rapid mode concentrates at wall', xlabel='Radius / wall radius', ylabel='Axial mean square / component maximum')
    axes[0, 0].legend(fontsize=8)
    growth = [wall['physical_mode'][0]['eigenvalue'], *wall['restricted_maximum_real'].values()]
    axes[0, 1].bar(['Original', 'Normal flow', 'Toroidal B', 'Both'], growth, color=['#bd6258', '#c79152', '#4687a6', '#448264'])
    axes[0, 1].set(yscale='log', title='Boundary restrictions remove rapid branch', xlabel='Imposed physical wall trace', ylabel='Largest real spatial growth rate (1/s)')
    axes[0, 1].text(.04, .04, 'Full constrained RHS: 9,203/s\nSlower growing modes remain', transform=axes[0, 1].transAxes, fontsize=8)
    for key, label in [('density', 'Density'), ('electron_pressure', 'Electron pressure'), ('electric', 'Full electric field')]:
        axes[0, 2].loglog([16, 32, 64], [a['order20']['norms']['final'][key]['relative_L2'] for _, a, _ in contacts], 'o-', label=label)
    axes[0, 2].set(title='Complete 100 microsecond contact transits', xlabel='Axial fluid points', ylabel='Relative physical signal error')
    axes[0, 2].legend(fontsize=8)
    for kind, label in [('wave', 'Wave B'), ('contact', 'Contact density')]:
        data = time[kind]; measure = data['physical_time_differences']['20']
        axes[1, 0].loglog(data['steps'][:-1], measure['differences'], 'o-', label=f"{label}; ratio {measure['ratios'][-1]:.2f}")
    axes[1, 0].set(title='Gauss time refinement, independent fields', xlabel='Steps in coarser full-period run', ylabel='Successive difference / physical signal norm')
    axes[1, 0].legend(fontsize=8)
    subset = [c for c in mode['cases'] if c['order'] == 20 and c['epsilon'] == 1e-5 and c['elements'] >= 4]
    for key, label in [('electric', 'Electric'), ('radial_velocity', 'Radial velocity'), ('azimuthal_B', 'Toroidal B rate'), ('current_rate', 'Current rate')]:
        axes[1, 1].loglog([2*c['elements'] for c in subset], [c['rate_errors'][key]['relative_L2'] for c in subset], 'o-', label=label)
    axes[1, 1].set(title='Independent growing-mode rate comparison', xlabel='Axial field cells', ylabel='Relative original-continuum rate error')
    axes[1, 1].legend(fontsize=8)
    data = time['contact']
    axes[1, 2].loglog(data['steps'], np.abs(data['continuous_energy_relative']), 'o-', label='Independent physical quadrature')
    axes[1, 2].loglog(data['steps'], np.abs(data['grid_energy_relative']), 'x--', label='Grid energy')
    axes[1, 2].set(title='Contact energy drift decreases with time step', xlabel='Steps per full transit', ylabel='|Total energy change| / initial total energy')
    axes[1, 2].legend(fontsize=8)
    for axis in axes.ravel(): axis.grid(alpha=.2)
    fig.suptitle('Wall correction enables complete contacts; finer mode and local-budget checks remain', fontsize=16, weight='bold')
    fig.savefig(destination/'convergence.png', dpi=150); fig.savefig(destination/'convergence.pdf'); plt.close(fig)
    if provenance(['examples.strong_wall_report']) != identity: raise RuntimeError('Report source changed')
    summary = dict(terminal=True, source_sha256=identity['source_sha256'], inputs=inputs, temporal=time,
        contact_mesh=[dict(elements=r['spec']['elements'], order16=a['order16'], order20=a['order20']) for r, a, _ in contacts],
        wave_mesh=[dict(elements=r['spec']['elements'], order16=a['order16'], order20=a['order20']) for r, a, _ in waves],
        limitations='Bounded periodic verification only. Growing-mode magnetic/current rate errors remain substantial at the finest tested mesh. Physical perturbation trajectories, local charge/enthalpy/Poynting interfaces, exact vacuum and formation remain open.')
    (destination/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    import sys
    run(sys.argv[1])
