"""Show rejected pressure-only states without concealing accepted/rejected times."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def figure(directory, failures, records):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    ax = axes[0]
    for record in failures:
        h = record['history_columns']; rejected = record['rejected_trial']; n = record['configuration']['steps']
        line, = ax.plot(np.array(h['time_s'])*1e6, h['min_ion_pressure_pa'], label=f'{n} planned steps')
        ax.plot(rejected['time_s']*1e6, rejected['min_ion_pressure_pa'], 'x', color=line.get_color(), markersize=8)
    ax.axhline(0, color='black', linewidth=.8)
    ax.set(title='nr128: no accepted negative ion pressure', xlabel='Time [µs]; × marks rejected trial', ylabel='Minimum ion pressure [Pa]')
    ax.grid(alpha=.18); ax.legend(fontsize=8)
    ax = axes[1]; steps = [80, 160, 320, 640]
    names = ['pressure_evolve_64']+[f'pressure_fine_time_{n}' for n in steps[1:]]
    errors = [records[n]['metrics']['ion_pressure_max_change_pa'] for n in names]
    ax.semilogx(steps, errors, 'o-', color='#b03a2e')
    ax.set_xticks(steps, labels=[str(n) for n in steps]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set(title='nr64: smaller time steps retain local error', xlabel='Accepted steps to 1 µs', ylabel='Maximum ion-pressure change [Pa]')
    ax.grid(alpha=.18)
    fig.suptitle('Pressure-only splitting limitation — conserved global energy does not ensure local accuracy')
    path = directory/'pressure_split_rejections.png'; fig.savefig(path, dpi=175); plt.close(fig)
    return path
