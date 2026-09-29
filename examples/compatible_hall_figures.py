"""Scientific figures from audited compatible Hall records."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def decorate(axes):
    for ax in np.asarray(axes).flat:
        ax.grid(True, which='both', alpha=.22)
        ax.legend(fontsize=8)


def figures(directory, records, comparison, legacy):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    paths = []
    nr = np.array([8, 16, 32, 64, 128])
    metrics = [records[f'spatial_n{n}']['metrics'] for n in nr]
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.1), constrained_layout=True)
    ax[0].loglog(nr, [r['error'] for r in legacy['first_annulus']], 's--', color='#ac4848', label='Earlier operator (f28714c)')
    ax[0].loglog(nr, [r['axis_g_rate_relative_L2'] for r in metrics], 'o-', label='Compatible curl / energy')
    ax[0].set(title='Regular-axis coefficient error', xlabel='Radial cells', ylabel='Relative first-annulus ġ error')
    for key, label in [('current_relative_L2', 'Current J'), ('induction_relative_L2', 'Physical ∂B/∂t'), ('g_rate_relative_L2', 'Cell ġ')]:
        ax[1].loglog(nr, [r[key] for r in metrics], 'o-', label=label)
    ax[1].set(title='Independent nonlinear snapshot', xlabel='Radial cells', ylabel='Relative volume L2 error')
    ax[2].loglog(nr, [r['local_energy_defect_relative_L2'] for r in metrics], 'o-', label='Local magnetic / Poynting difference')
    ax[2].loglog(nr, metrics[0]['local_energy_defect_relative_L2']*(8/nr)**2, '--', color='gray', label='Second-order guide')
    ax[2].set(title='Local Hall transfer defect', xlabel='Radial cells', ylabel='Relative volume L2 difference')
    decorate(ax)
    fig.suptitle('Compatible Hall spatial verification • finite-amplitude variable-density fields')
    path = directory/'spatial.png'; fig.savefig(path, dpi=160); plt.close(fig); paths.append(path)

    nr = [8, 16, 32]
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.1), constrained_layout=True)
    for branch, label, color in [(0, 'Slow: 0.5 rad', '#b55d22'), (1, 'Intermediate: 2π', '#308558'), (2, 'Fast: 2π', '#376bab')]:
        data = [records[f'wave_b{branch}_n{n}'] for n in nr]
        ax[0].loglog(nr, [r['errors']['magnetic_relative_L2'] for r in data], 'o-', color=color, label=label)
        ax[0].loglog(nr, [r['initial_errors']['magnetic_relative_L2'] for r in data], ':', color=color)
        ax[1].loglog(nr, [r['errors']['velocity_relative_L2'] for r in data], 'o-', color=color, label=label)
        ax[2].loglog(nr, [r['errors']['g_cell_relative_L2'] for r in data], 'o-', color=color, label=label)
    ax[0].plot([], [], ':', color='gray', label='Dotted: initial representation error')
    ax[0].set(title='Physical magnetic perturbation', xlabel='Radial cells', ylabel='Relative volume L2 error')
    ax[1].set(title='Physical velocity', xlabel='Radial cells', ylabel='Relative volume L2 error')
    ax[2].set(title='Toroidal coefficient cell averages', xlabel='Radial cells', ylabel='Relative volume L2 error')
    decorate(ax)
    fig.suptitle('Actual compressible Hall evolution • linear continuum references • no bias in error denominator')
    path = directory/'waves.png'; fig.savefig(path, dpi=160); plt.close(fig); paths.append(path)

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.1), constrained_layout=True)
    for family, counts in [('wave', (8, 16, 32)), ('nonlinear', (4, 8, 16))]:
        values = [comparison['fixed_grid_time'][f'{family}_time_{n}'] for n in counts]
        ax[0].loglog(counts, [v['magnetic_relative_to_reference_change'] for v in values], 'o-', label=family+' magnetic')
        ax[0].loglog(counts, [v['momentum_relative_to_reference_change'] for v in values], 's--', label=family+' momentum')
    ax[0].set(title='Fixed-grid time refinement', xlabel='Accepted steps', ylabel='Error / fine-reference change')
    values = [comparison['nonlinear_nested_mesh'][f'nonlinear_n{n}'] for n in nr]
    for key, label in [('magnetic', 'Magnetic'), ('momentum', 'Momentum'), ('density', 'Density')]:
        ax[1].loglog(nr, [v[key+'_relative_to_reference_change'] for v in values], 'o-', label=label)
    ax[1].set(title='Nonlinear evolution vs restricted nr64', xlabel='Radial cells', ylabel='Difference / reference change')
    nr = [8, 16, 32, 64]
    values = [records[f'decay_n{n}_s64']['metrics'] for n in nr]
    for key, label in [('g_relative_L2', 'Toroidal field coefficient'), ('axis_g_relative_error', 'First annulus'), ('physical_heat_relative_error', 'Physical Ohmic heat')]:
        ax[2].loglog(nr, [abs(v[key]) for v in values], 'o-', label=label)
    ax[2].set(title='Independent conducting Bessel decay', xlabel='Radial cells', ylabel='Relative error')
    decorate(ax)
    fig.suptitle('Separate time, nonlinear mesh and resistive controls • unchanged initial conditions')
    path = directory/'controls.png'; fig.savefig(path, dpi=160); plt.close(fig); paths.append(path)
    return paths
