"""Measured bounded coupled progress, including retained accuracy limitations."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from .regular_plasma_reference import initialize


def report(directory, destination):
    directory = Path(directory); destination = Path(destination); destination.mkdir(parents=True, exist_ok=True)
    records = [json.loads(p.read_text()) for p in sorted(directory.glob('case_*.json'))]
    audit = json.loads((directory/'audit.json').read_text())
    assert audit['terminal'] and len(audit['records']) == len(records)
    contact = [records[i] for i in (0, 1, 2)]
    wave = [records[i] for i in (3, 4, 5)]
    evolved = [records[i] for i in (12, 13, 14)]
    times = [records[i] for i in (15, 12, 16, 17)]
    b = [np.load(directory/r['arrays_file'])['final_magnetic'] for r in times]
    _, _, _, system, _ = initialize(times[0]['spec'])
    h = system.geometry.hessian
    differences = [float(np.sqrt((b[i]-b[i+1])@(h@(b[i]-b[i+1]))/(b[-1]@(h@b[-1])))) for i in range(3)]
    ratios = [differences[i]/differences[i+1] for i in range(2)]
    spectra = [records[i] for i in (19, 20)]
    summary = dict(source_sha256=audit['source_sha256'], cases=len(records),
        terminal_audit=True, wave_period=evolved[0]['duration'],
        contact_rate_finest=contact[-1]['independent'], wave_rate_finest=wave[-1]['independent'],
        full_wave_finest=evolved[-1]['final'], magnetic_time_differences=differences,
        magnetic_time_ratios=ratios, stationary_spectra=spectra,
        max_rate_residual=max(r['residual'] for r in records if r['kind']=='rate'),
        max_trajectory_residual=max(r['max_residual'] for r in records if r['kind']=='evolution'),
        energy_change_over_wave_energy=[r['energy_change_over_wave_energy'] for r in evolved],
        long_contact_trajectory_complete=False, full_residual_refinement_adopted=False,
        physical_interface_energy_budget_complete=False, formation_verified=False)
    (destination/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    for key, label in [('density_rate', 'density'), ('electron_pressure_rate', 'electron pressure'),
                       ('ion_pressure_rate', 'ion pressure'), ('electric', 'electric field')]:
        axes[0, 0].loglog([r['spec']['elements'] for r in contact],
                         [r['independent']['norms'][key]['relative_L2'] for r in contact], 'o-', label=label)
    for key, label in [('velocity_rate', 'velocity'), ('magnetic_rate', 'magnetic'), ('current_rate', 'current'), ('electric', 'electric')]:
        axes[0, 1].loglog([r['spec']['elements'] for r in wave],
                         [r['independent']['norms'][key]['relative_L2'] for r in wave], 'o-', label=label)
    for key, label in [('density', 'density / both pressures'), ('velocity', 'velocity'), ('magnetic', 'magnetic'), ('current', 'current')]:
        axes[0, 2].loglog([r['spec']['elements'] for r in evolved],
                         [r['final'][key]['relative_L2'] for r in evolved], 'o-', label=label)
    axes[1, 0].loglog([r['steps'] for r in times[:-1]], differences, 'o-', label='successive magnetic difference')
    axes[1, 0].text(.05, .08, 'halving-step ratios: '+', '.join(f'{v:.3f}' for v in ratios), transform=axes[1, 0].transAxes)
    for branch, indices in [(0, (8, 21, 22)), (1, (9, 23, 24))]:
        rows = [records[i] for i in indices]
        axes[1, 1].loglog([r['spec']['elements'] for r in rows],
            [r['independent']['norms']['magnetic_rate']['relative_L2'] for r in rows], 'o-', label=f'branch {branch}: magnetic rate')
    values = np.load(directory/records[20]['arrays_file'])['eigenvalues']
    axes[1, 2].plot(values.real, values.imag/1e6, '.', markersize=3)
    axes[1, 2].set_xlabel('real part (1/s)'); axes[1, 2].set_ylabel('imaginary part (Mrad/s)')
    axes[1, 2].text(.04, .03, f"energy symmetry defect {spectra[-1]['energy_symmetric_defect']:.2e}\nmax real / radius {spectra[-1]['maximum_real_over_radius']:.2e}", transform=axes[1, 2].transAxes)
    titles = ['Exact moving contact · instantaneous rates', 'Fast finite-mass wave · instantaneous rates',
              'One complete fast-wave period', 'Full coupled RK4 time refinement',
              'Slower branches · retained rate errors', 'Full separate stationary spectrum · 360 DOFs']
    for ax, title in zip(axes.flat, titles):
        ax.set_title(title); ax.grid(alpha=.2)
    for ax in (axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 1]):
        ax.set_xlabel('radial fluid elements'); ax.set_ylabel('relative physical L2 error'); ax.legend(fontsize=8)
    axes[1, 0].set_xlabel('steps per wave period'); axes[1, 0].set_ylabel('relative magnetic difference')
    fig.suptitle('Coupled moving plasma: bounded verification progress\nFilled periodic cylinder · physical masses · no added heat or diffusion', fontsize=17)
    fig.supxlabel('Long contact evolution, strict residual robustness, physical interface energy budgets and exact-vacuum formation remain open.', fontsize=10)
    fig.savefig(destination/'convergence.png', dpi=165)
    fig.savefig(destination/'convergence.pdf')
    plt.close(fig)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('directory'); parser.add_argument('destination')
    args = parser.parse_args(); print(json.dumps(report(args.directory, args.destination), indent=2))
