"""Audit complete Hall spatial records and plot accuracy separately from work."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.toroidal_state import radius_squared
from .hall_work_cases import protocol
from .hall_work_measurements import provenance
from .verification_records import write_record


def run(directory='outputs/hall_work_verification'):
    directory = Path(directory)
    data = json.loads((directory/'measurements.json').read_text())
    if not data['completed'] or data['provenance'] != provenance():
        raise ValueError('Incomplete or stale Hall spatial study')
    expected = protocol()
    if [r['configuration'] for r in data['cases']] != expected:
        raise ValueError('Hall spatial case coverage/configuration mismatch')
    by_name, replay_errors = {}, []
    for row in data['cases']:
        case = row['configuration']
        name = case['name']
        if not row['completed'] or row['provenance'] != data['provenance'] or row['kind'] != 'instantaneous_spatial_audit':
            raise ValueError('Partial or misclassified Hall record')
        if json.loads((directory/(name+'.json')).read_text()) != row:
            raise ValueError('Collected/individual Hall record disagreement')
        path = directory/row['fields']['file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['fields']['sha256']:
            raise ValueError('Hall lossless archive identity mismatch')
        with np.load(path) as fields:
            if not all(np.isfinite(fields[key]).all() for key in fields.files):
                raise ValueError('Nonfinite Hall archive values')
            mesh = MHDMesh(case['reference']['radius'], case['reference']['length']/2, case['nr'], case['nz'])
            if fields['g'].shape != (mesh.nr, mesh.nz) or fields['psi'].shape != (mesh.nr+1, mesh.nz+1):
                raise ValueError('Hall archive geometry mismatch')
            br, bz = mesh.cell_fields(fields['psi'])
            dbr, dbz = mesh.cell_fields(fields['psi_rate'])
            rate = (br*dbr+bz*dbz+radius_squared(mesh)*fields['g']*fields['g_rate'])/MU0
            # The uniform bias makes each physical B dot Bdot term O(epsilon),
            # while Hall work is O(epsilon²). Use the actual summation condition
            # scale for this independent floating-point replay, not the latter.
            scale = np.sum(mesh.volumes*(abs(br*dbr)+abs(bz*dbz)+abs(radius_squared(mesh)*fields['g']*fields['g_rate'])))/MU0
            difference = abs(np.sum(mesh.volumes*rate)-row['metrics']['global_magnetic_work_w'])
            replay_errors.append(float(difference/scale) if scale > 0 else 0.)
            if difference > 64*np.finfo(float).eps*scale:
                raise ValueError('Archived physical magnetic work disagrees with report')
            if not np.array_equal(fields['psi_rate'][[0, -1]], np.zeros((2, mesh.nz+1))) or not np.array_equal(fields['radial_poynting'][[0, -1]], np.zeros((2, mesh.nz))):
                raise ValueError('Conducting wall / regular axis trace mismatch')
            if not np.array_equal(fields['axial_poynting'][:, 0], fields['axial_poynting'][:, -1]) or np.any(fields['density_quadrature'] <= 0):
                raise ValueError('Periodic energy or positive-carrier assumption violated')
        by_name[name] = row
    if len(list(directory.glob('*.npz'))) != len(expected):
        raise ValueError('Unexpected untracked Hall case archive')
    base = [by_name[f'nonlinear_n{n}'] for n in (8, 16, 32, 64, 128)]
    nr = np.array([r['configuration']['nr'] for r in base])
    values = lambda key: np.array([r['metrics'][key] for r in base])
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    ax = axes[0, 0]
    for key, label in (('induction_relative_L2', 'Magnetic rate at physical points'), ('current_relative_L2', 'Compatible current'), ('psi_rate_max_relative', 'Corner flux rate, max norm')):
        ax.loglog(nr, values(key), 'o-', label=label)
    ax.set(title='Global spatial error decreases', ylabel='Relative error')
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    ax.semilogx(nr, 100*values('axis_g_rate_relative_L2'), 'o-', color='#b54736')
    ax.set(title='Axis consistency gate fails', ylabel='First-annulus g rate error (%)', ylim=(0, 25))
    ax.text(.97, .12, 'Not enabled in the fluid solver', transform=ax.transAxes, ha='right', fontsize=9, color='#b54736')
    ax = axes[1, 0]
    ax.loglog(nr, values('local_energy_defect_relative_L2'), 'o-', color='#7556a6')
    ax.set(title='Local energy transfer remains approximate', ylabel='Magnetic rate minus flux rate, relative L2')
    ax = axes[1, 1]
    ax.semilogx(nr, values('global_magnetic_work_relative')*1e18, 'o-', color='#237a57')
    ax.set(title='Global Hall work cancels', ylabel='Residual / integral |J| |E|  (×10⁻¹⁸)')
    for ax in axes.flat:
        ax.set_xlabel('Radial cells (axial cells = 2 × radial cells)')
        ax.grid(alpha=.22)
        ax.set_xticks(nr, labels=[str(n) for n in nr])
    fig.suptitle('Experimental nonlinear Hall operator: conservation does not certify the axis', fontsize=13)
    first = directory/'hall_spatial_accuracy.png'
    fig.savefig(first, dpi=170); plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for n, style in ((16, '--'), (128, '-')):
        with np.load(directory/f'nonlinear_n{n}.npz') as f:
            z = .5*(f['axial_faces'][1:]+f['axial_faces'][:-1])
            axes[0].plot(z, f['g_rate'][0], style, label=f'Computed, nr={n}')
            if n == 128:
                axes[0].plot(z, f['reference_g_rate'][0], ':', color='black', label='Continuum annular mean')
    axes[0].set(title='First-annulus waveform', xlabel='z (m)', ylabel='g rate (T m⁻¹ s⁻¹)')
    axes[0].legend(fontsize=8)
    for branch in (0, 1, 2):
        for sign, style in ((1, '-'), (-1, ':')):
            y = [by_name[f'wave_b{branch}_k{sign}_n{n}']['metrics']['induction_relative_L2'] for n in (8, 16, 32)]
            axes[1].loglog([8, 16, 32], y, 'o'+style, label=f'Branch {branch}, k {sign:+d}')
    axes[1].set(title='Finite-pressure Hall induction only', xlabel='Radial cells', ylabel='Relative physical-volume rate error')
    axes[1].legend(fontsize=7, ncol=2)
    names = ['wave_amplitude_0.001', 'wave_b2_k1_n32', 'wave_amplitude_0.1']
    amplitude = [by_name[name]['configuration']['reference']['amplitude'] for name in names]
    quadratic = [by_name[name]['metrics']['quadratic_induction_rms_t_s'] for name in names]
    symmetric = [by_name[name]['metrics']['quadratic_symmetric_difference_rms_t_s'] for name in names]
    axes[2].loglog(amplitude, quadratic, 'o-', color='#7556a6', label='Direct quadratic evaluation')
    axes[2].loglog(amplitude, symmetric, 'x--', color='#b54736', label='Full-field cancellation')
    axes[2].legend(fontsize=8)
    axes[2].set(title='Resolved quadratic Hall contribution', xlabel='Mode amplitude (m/s)', ylabel='Quadratic magnetic rate RMS (T/s)')
    for ax in axes:
        ax.grid(alpha=.22)
    second = directory/'hall_axis_and_modes.png'
    fig.savefig(second, dpi=170); plt.close(fig)
    q3 = by_name['nonlinear_n32']['metrics']
    q5 = by_name['quadrature_q5_n32']['metrics']
    q7 = by_name['quadrature_q7_n32']['metrics']
    summary = dict(completed=True, completed_spatial_cases=len(expected), time_evolution_verified=False,
                   enabled_in_fluid_solver=False, source_sha256=data['provenance']['source_sha256'],
                   baseline={key: values(key).tolist() for key in ('current_relative_L2', 'induction_relative_L2', 'g_rate_cell_relative_L2', 'axis_g_rate_relative_L2', 'psi_rate_max_relative', 'local_energy_defect_relative_L2', 'global_magnetic_work_relative')},
                   max_global_work_relative=max(row['metrics']['global_magnetic_work_relative'] for row in data['cases']),
                   max_archive_work_replay_backward_error=max(replay_errors),
                   max_reference_quadrature_change=max(row['metrics']['cell_reference_quadrature_relative_change'] for row in data['cases']),
                   induction_error_change_q3_q5=abs(q3['induction_relative_L2']-q5['induction_relative_L2']),
                   induction_error_change_q5_q7=abs(q5['induction_relative_L2']-q7['induction_relative_L2']),
                   quadratic_rate_amplitudes=amplitude, quadratic_rates_t_s=quadratic,
                   quadratic_symmetric_difference_rates_t_s=symmetric,
                   quadratic_symmetric_difference_errors=[by_name[name]['metrics']['quadratic_symmetric_difference_relative_error'] for name in names],
                   quadratic_scaling_powers=(np.diff(np.log(quadratic))/np.diff(np.log(amplitude))).tolist(),
                   report_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   figures={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (first, second)})
    write_record(directory/'metrics.json', summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    run()
