"""Static scientific figures of measured warm Hall convergence and limitations."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def finish(fig, path):
    fig.savefig(path, dpi=175); plt.close(fig)
    return path


def convergence(directory, records):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    grids = np.array([4, 8, 16, 32])
    for branch in range(3):
        ax = axes[0, branch]
        for key, label in (('electron_pressure_relative_perturbation_L2', 'Electron pressure'),
                            ('ion_pressure_relative_perturbation_L2', 'Ion pressure'), ('magnetic_relative_L2', 'Magnetic field')):
            ax.loglog(grids, [records[f'wave_b{branch}_{n}']['metrics'][key] for n in grids], 'o-', label=label)
        ax.set(title=f'Warm wave branch {branch}: phase 0.5 rad', xlabel='Radial cells (axial = 2×)', ylabel='Physical relative L2 perturbation error')
        ax.set_xticks(grids, labels=grids.astype(str)); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(which='both', alpha=.18); ax.legend(fontsize=8)
    grids = np.array([8, 16, 32, 64, 128]); ax = axes[1, 0]
    for key, label in (('energy_error_relative_to_pressure_work', 'Conservative energy rate'),
                       ('local_ion_defect_relative_to_pressure_work', 'Local ion-energy defect'),
                       ('axis_energy_error_relative_to_pressure_work', 'First-annulus energy rate')):
        ax.loglog(grids, [records[f'pressure_spatial_{n}']['metrics'][key] for n in grids], 'o-', label=label)
    ax.set(title='Independent nonlinear pressure work', xlabel='Radial cells', ylabel='L2 error / pressure-work L2')
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    values = [records[f'trace_diffusion_{n}']['metrics'] for n in grids]
    ax.loglog(grids, [abs(v['axis_instantaneous_g_rate_relative_signed_error']) for v in values], 'o-', label='Instantaneous axis derivative')
    ax.loglog(grids, [v['evolved_g_relative_L2'] for v in values], 'o-', label='Evolved g (volume norm)')
    ax.loglog(grids, [v['evolved_magnetic_coefficient_relative_L2'] for v in values], 'o-', label='Evolved g (magnetic weight)')
    ax.set(title='Standalone resistive boundary test', xlabel='Radial cells', ylabel='Relative error; decay exponent 0.02')
    ax.legend(fontsize=8)
    ax = axes[1, 2]; ns = grids[:-1]
    ax.loglog(ns, [records[f'pressure_evolve_{n}']['metrics']['ion_pressure_rms_change_pa'] for n in ns], 'o-', label='Volume RMS')
    ax.loglog(ns, [records[f'pressure_evolve_{n}']['metrics']['ion_pressure_max_change_pa'] for n in ns], 'o-', label='Maximum local error')
    ax.set(title='Stationary-ion pressure subproblem', xlabel='Radial cells', ylabel='Spurious ion-pressure change [Pa] at 1 µs')
    ax.legend(fontsize=8)
    for ax in axes[1]:
        ax.grid(which='both', alpha=.18); ax.set_xticks(grids, labels=grids.astype(str))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.suptitle('Warm Hall verification — bounded periodic cylinder; no formation claim')
    return finish(fig, directory/'warm_hall_convergence.png')


def controls(directory, records, arrays, comparisons):
    from .warm_hall_reference import WarmHallWave
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8), constrained_layout=True)
    for ax, branch in zip(axes[0], (0, 2)):
        name = f'wave_b{branch}_32'; record = records[name]; a = arrays[name]
        ref = WarmHallWave(**record['configuration']['reference'])
        i = 8; radius = a['radial_faces'][i]
        z = a['axial_faces'][:-1]; smooth = np.linspace(z[0], a['axial_faces'][-1], 600)
        true = ref.fields(radius, smooth, record['actual_time_s'])['electron_pressure']
        actual = a['final_electron_pressure'][i]-ref.electron_pressure
        ax.plot(smooth, true, color='#333333', label='Independent linear wave')
        ax.plot(z, actual, '.', color='#a144c7', label='Electron pressure nodes')
        ax.set(title=f'Branch {branch}, nr=32, r/R=0.25', xlabel='z [m]', ylabel='Electron-pressure perturbation [Pa]')
        ax.grid(alpha=.18); ax.legend(fontsize=8)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax = axes[1, 0]
    for label, value in comparisons['time'].items():
        steps = [records[v['coarse']]['accepted_steps'] for v in value['comparisons']]
        error = [v['electron_pressure']['rms_difference_pa'] for v in value['comparisons']]
        ax.loglog(steps, error, 'o-', label=label)
    ax.set(title='Fixed-grid temporal differences', xlabel='Accepted steps', ylabel='Q1 electron-pressure difference [Pa]')
    ax.grid(which='both', alpha=.18); ax.legend(fontsize=8)
    ax = axes[1, 1]
    for n in (8, 16, 32):
        h = records[f'nonlinear_{n}']['history_columns']; initial = h['total_energy_j'][0]
        ax.plot(np.array(h['time_s'])*1e6, (np.array(h['total_energy_j'])-initial)/initial, label=f'nr={n}')
    ax.set(title='Nonlinear coupled total-energy budget', xlabel='Time [µs]', ylabel='Relative total-energy change')
    ax.grid(alpha=.18); ax.legend(fontsize=8)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    fig.suptitle('Species accuracy and conservative accounting are separate checks')
    return finish(fig, directory/'warm_hall_controls.png')
