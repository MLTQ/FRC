"""Validate complete Hall records and render measured scientific accuracy figures."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
from scipy.constants import elementary_charge as E, proton_mass as MP, mu_0 as MU0
from .hall_wave_cases import protocol
from .hall_wave_measurements import provenance
from .hall_wave_reference import frequency
from .verification_records import write_record


def validate(directory):
    data = json.loads((directory/'measurements.json').read_text())
    parallel = json.loads((directory/'parallel.json').read_text())
    if not data['completed'] or data['provenance'] != provenance():
        raise ValueError('Incomplete or stale Hall measurement sources')
    cases = protocol()
    if [c.name for c in cases] != [r['name'] for r in data['cases']]:
        raise ValueError('Hall case coverage differs from declared protocol')
    for case, row in zip(cases, data['cases']):
        if not row['completed'] or row['configuration'] != case.configuration() or row['provenance'] != data['provenance']:
            raise ValueError('Hall case identity differs')
        if json.loads((directory/(case.name+'.json')).read_text()) != row:
            raise ValueError('Collected and individual Hall records differ')
        path = directory/row['fields']['file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['fields']['sha256']:
            raise ValueError('Lossless Hall field archive differs')
        with np.load(path, allow_pickle=False) as fields:
            if any(not np.isfinite(fields[k]).all() for k in fields.files):
                raise ValueError('Nonfinite archived Hall field')
        if len(row['history']) != case.steps or any(h['heat_j'] < 0 for h in row['history']):
            raise ValueError('Incomplete history or negative physical heat')
    expected = [(m, kd, h) for m in ('hall', 'two_fluid') for kd in (.1, 1., 10., 1000.) for h in (-1, 1)]
    if not parallel['completed'] or [(r['model'], r['k_di'], r['helicity']) for r in parallel['cases']] != expected:
        raise ValueError('Incomplete parallel-wave protocol')
    for name, sha in parallel['source_files'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha:
            raise ValueError('Stale parallel source')
    return data, parallel


def plots(directory, data, parallel):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    fig.suptitle('Linear Hall verification: cylinder geometry, both helicities and finite electron mass', fontsize=15)
    fig.supxlabel('Uniform n = 10¹⁸ m⁻³, B₀ = 1 mT; cylinder R = 1 m for Hall, 0.01 m for finite mass; L = 2R', fontsize=9)
    colors = {-1: '#006c8c', 1: '#bf4b1c'}
    for model, linestyle in (('hall', '--'), ('two_fluid', '-')):
        for h in (-1, 1):
            rows = [r for r in data['cases'] if r['configuration']['control'] == 'base' and r['configuration']['model'] == model and r['configuration']['helicity'] == h]
            nr = [r['configuration']['nr'] for r in rows]
            label = f'{model.replace("_", " ")}, h={h:+d}'
            axes[0, 0].loglog(nr, [r['final_errors']['magnetic_relative_volume_L2'] for r in rows], linestyle, marker='o', color=colors[h], label=label)
            axes[0, 1].loglog(nr, [abs(r['phase_frequency_relative_error']) for r in rows], linestyle, marker='o', color=colors[h], label=label)
    axes[0, 0].set(title='Actual cylinder field error', xlabel='Radial cells', ylabel='Relative volume L2 error')
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set(title='Measured cylinder phase error', xlabel='Radial cells', ylabel='Absolute relative frequency error')
    for h in (-1, 1):
        rows = [r for r in data['cases'] if r['configuration']['model'] == 'two_fluid' and r['configuration']['nr'] == 32 and r['configuration']['helicity'] == h and (r['configuration']['control'] == 'base' or r['configuration']['control'].startswith('time'))]
        rows.sort(key=lambda r: r['configuration']['steps'])
        axes[0, 2].loglog([r['configuration']['steps'] for r in rows], [r['time_relative_energy_norm_error'] for r in rows], 'o-', color=colors[h], label=f'h={h:+d}')
    axes[0, 2].set(title='Time error against spatial exponential', xlabel='Steps over two radians', ylabel='Relative complete-energy norm')
    axes[0, 2].legend(fontsize=9)
    kd = np.geomspace(.05, 2000, 300)
    di = np.sqrt(MP/(MU0*1e18*E**2))
    oi = E*.001/MP
    for model, linestyle in (('hall', '--'), ('two_fluid', '-')):
        for h in (-1, 1):
            values = [frequency(1e18, MP, .001, 3e-5, value/di, 0., h, model).real/oi for value in kd]
            axes[1, 0].loglog(kd, values, linestyle, color=colors[h])
            rows = [r for r in parallel['cases'] if r['model'] == model and r['helicity'] == h]
            axes[1, 0].loglog([r['k_di'] for r in rows], [r['frequency_real_s-1']*(1+r['refinements'][-1]['phase_frequency_relative_error'])/oi for r in rows], 'o', color=colors[h], markersize=4)
    axes[1, 0].set(title='Parallel dispersion: lines exact, dots evolved', xlabel='k di', ylabel='Frequency / ion gyrofrequency')
    rows = [r for r in data['cases'] if r['configuration']['model'] == 'two_fluid' and r['configuration']['nr'] == 32 and r['configuration']['control'] == 'base']
    for row in rows:
        h = row['configuration']['helicity']
        hist = row['history']
        axes[1, 1].plot([v['time_s']/row['end_time_s'] for v in hist], [v['relative_energy_balance'] for v in hist], color=colors[h], label=f'h={h:+d}')
    axes[1, 1].set(title='Energy + heat balance (not phase accuracy)', xlabel='Time / final time', ylabel='Residual / initial energy')
    axes[1, 1].legend(fontsize=9)
    for h in (-1, 1):
        row = next(r for r in parallel['cases'] if r['model'] == 'two_fluid' and r['k_di'] == 1000. and r['helicity'] == h)
        axes[1, 2].loglog([s['steps'] for s in row['refinements']], [s['relative_state_error'] for s in row['refinements']], 'o-', color=colors[h], label=f'h={h:+d}')
    axes[1, 2].set(title='Independent parallel-wave time error', xlabel='Steps over two radians', ylabel='Relative full-state error')
    axes[1, 2].legend(fontsize=9)
    for ax, ticks in ((axes[0, 0], [8, 16, 32, 64]), (axes[0, 1], [8, 16, 32, 64]),
                      (axes[0, 2], [128, 512, 2048, 8192]), (axes[1, 2], [40, 80, 160])):
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    for ax in axes.flat:
        ax.grid(alpha=.18)
    fig.savefig(directory/'hall_accuracy.png', dpi=160)
    plt.close(fig)
    case = next(c for c in protocol() if c.name == 'two_fluid_h+1_nr64_base')
    solver, reference, _ = case.build()
    row = next(r for r in data['cases'] if r['name'] == case.name)
    fields = np.load(directory/row['fields']['file'], allow_pickle=False)
    r, z = solver.mesh.r[:, None], np.linspace(0., case.length, 201)[None, :]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout='constrained')
    fig.suptitle('Regular cylindrical electron-branch wave — linear reference, no formation claim')
    limits = {q: np.max(abs(solver.sample(fields['initial_modes'], r, z, q)[..., 1])) for q in ('magnetic', 'velocity')}
    for column, (key, time, title) in enumerate((('initial_modes', 0., 'Initial'), ('final_modes', row['end_time_s'], 'Evolved'))):
        for axis, quantity, label, scale in ((axes[0, column], 'magnetic', 'B theta (µT)', 1e6), (axes[1, column], 'velocity', 'v theta (m/s)', 1.)):
            values = solver.sample(fields[key], r, z, quantity)[..., 1]*scale
            limit = limits[quantity]*scale
            im = axis.pcolormesh(z.ravel()*1e3, r.ravel()*1e3, values, cmap='RdBu_r', vmin=-limit, vmax=limit, shading='auto', rasterized=True)
            axis.set(title=title+' '+label, xlabel='z (mm)', ylabel='r (mm)')
            fig.colorbar(im, ax=axis, pad=.02)
    fig.savefig(directory/'hall_wave.png', dpi=160)
    plt.close(fig)


def run(directory='outputs/hall_verification'):
    directory = Path(directory)
    data, parallel = validate(directory)
    plots(directory, data, parallel)
    metrics = {'completed': True, 'cylinder_case_count': len(data['cases']), 'parallel_mode_count': len(parallel['cases']),
               'source_sha256': data['provenance']['source_sha256'],
               'maximum_energy_balance_relative': max(abs(h['relative_energy_balance']) for r in data['cases'] for h in r['history']),
               'maximum_equation_relative_residual': max(h['equation_relative_residual'] for r in data['cases'] for h in r['history']),
               'maximum_physical_energy_relative_error': max(r['physical_vs_modal_relative_error'] for r in data['cases']),
               'maximum_error_quadrature_absolute_change': max(r['error_quadrature_absolute_change'] for r in data['cases']),
               'parallel_time_ratio_range': [min(min(r['time_error_ratios']) for r in parallel['cases']), max(max(r['time_error_ratios']) for r in parallel['cases'])],
               'report_source_files': {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in ('examples/verify_hall.py', 'examples/verification_records.py')},
               'figures': {name: hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in ('hall_accuracy.png', 'hall_wave.png')}}
    write_record(directory/'metrics.json', metrics)
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    run()
