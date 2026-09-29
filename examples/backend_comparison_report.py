"""Replay external-engine archives and plot the bounded convergence comparison."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import c, epsilon_0, mu_0, proton_mass
from .backend_wave_reference import ParallelWave, relative_error, vacuum_field


def run(directory='outputs/engine_comparison'):
    directory = Path(directory)
    rows = []
    paths = sorted((directory/'runs').glob('*/measurement.json'))
    expected = {f'hall_h{h}_n{n}' for h in (-1, 1) for n in (32, 64, 128)}
    expected |= {f'vacuum_n{n}' for n in (32, 64, 128)}
    expected |= {f'hall_h{h}_n128_half_dt' for h in (-1, 1)}
    expected |= {f'hall_h{h}_n128_quarter_dt' for h in (-1, 1)}
    expected |= {'hall_h1_n128_ppc64', 'hall_h1_n128_small_amplitude',
                 'hall_h1_n128_floor_control', 'vacuum_n128_half_dt'}
    if {p.parent.name for p in paths} != expected:
        raise RuntimeError('The declared 17-case protocol is incomplete or contains extra cases')
    for path in paths:
        row = json.loads(path.read_text())
        archive = path.parent/'wave_samples.npz'
        assert row['completed']
        assert row['actual_steps'] == row['steps']
        np.testing.assert_allclose(row['actual_time_s'], row['stop_s'], rtol=1e-12, atol=0.)
        assert hashlib.sha256(archive.read_bytes()).hexdigest() == row['archive_sha256']
        for name, digest in row['source_files'].items():
            assert hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() == digest
        with np.load(archive, allow_pickle=False) as data:
            if row['kind'] == 'hall':
                ref = ParallelWave(**row['parameters'])
                np.testing.assert_allclose(data['Ex_initial'], ref.omega/ref.k
                    *ref.magnetic(data['Ex_z'], 0.)[1], rtol=1e-10, atol=1e-18)
                np.testing.assert_allclose(data['Ey_initial'], -ref.omega/ref.k
                    *ref.magnetic(data['Ey_z'], 0.)[0], rtol=1e-10, atol=1e-18)
                actual = np.array([data['Bx_final'], data['By_final']])
                error = relative_error(actual, ref.magnetic(data['Bx_z'], row['stop_s']))
                np.testing.assert_allclose(error, row['relative_magnetic_error'], rtol=1e-8, atol=1e-12)
                energies = []
                for endpoint in ('initial', 'final'):
                    np.testing.assert_allclose(data['Bz_'+endpoint], ref.bias, rtol=0., atol=1e-14)
                    u2 = sum(data['particle_u'+axis+'_'+endpoint]**2 for axis in 'xyz')
                    weight = data['particle_weight_'+endpoint]
                    # Independent physical count and cold transverse-energy replay.
                    np.testing.assert_allclose(weight.sum(), ref.density*ref.length, rtol=2e-12)
                    e_ion = (weight*proton_mass*u2/(np.sqrt(1+u2/c**2)+1)).sum()
                    e_b = sum((data['B'+axis+'_'+endpoint]**2).mean()
                              for axis in 'xy')*ref.length/(2*mu_0)
                    energies.append(e_ion+e_b)
                np.testing.assert_allclose(energies, row['perturbation_energy_J_m-2'], rtol=1e-8)
                energy_change = (energies[-1]-energies[0])/energies[0]
                row['audited_particle_number_relative_change'] = float(
                    data['particle_weight_final'].sum()/data['particle_weight_initial'].sum()-1)
            else:
                for component in ('By', 'Ex'):
                    exact = vacuum_field(data[component+'_z'], row['stop_s'],
                        amplitude=row['parameters']['amplitude_T'], electric=component == 'Ex')
                    error = relative_error(data[component+'_final'], exact)
                    np.testing.assert_allclose(error, row[component+'_relative_error'], rtol=1e-8, atol=1e-12)
                energies = [sum(np.mean(data['B'+a+'_'+s]**2)/(2*mu_0)
                    +epsilon_0*np.mean(data['E'+a+'_'+s]**2)/2 for a in 'xyz')
                    for s in ('initial', 'final')]
                energy_change = (energies[-1]-energies[0])/energies[0]
            np.testing.assert_allclose(energy_change, row['relative_endpoint_energy_change'], rtol=1e-7, atol=1e-12)
        history = np.loadtxt(path.parent/'field_energy.txt')
        assert len(history) == row['steps']+1
        np.testing.assert_array_equal(history[:, 0], np.arange(row['steps']+1))
        np.testing.assert_allclose(history[-1, 1], row['stop_s'], rtol=1e-12, atol=0.)
        if row['kind'] == 'hall':
            particle = np.loadtxt(path.parent/'particle_energy.txt')
            np.testing.assert_allclose(history[:, :2], particle[:, :2], rtol=1e-12, atol=0.)
            scale = row['perturbation_energy_J_m-2'][0]
            drift = ((history[:, 4]-history[0, 4])+(particle[:, 2]-particle[0, 2]))/scale
            # Native text has 14 digits after the decimal; B0 subtraction loses
            # relative precision compared with the independent endpoint archive.
            bound = 1e-14*(abs(history[:, 4]).max()+abs(history[0, 4])
                          +abs(particle[:, 2]).max()+abs(particle[0, 2]))/scale
        else:
            drift = (history[:, 2]-history[0, 2])/history[0, 2]
            bound = 3e-14
        assert abs(drift[-1]-row['relative_endpoint_energy_change']) <= bound+1e-10
        row['native_history_max_relative_energy_change'] = float(np.max(abs(drift)))
        row['native_history_print_precision_bound'] = float(bound)
        row['native_artifact_sha256'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.parent.iterdir() if p.is_file()}
        row['case'] = path.parent.name
        rows.append(row)
    convergence = {}
    for kind, helicity in [('hall', -1), ('hall', 1), ('vacuum', None)]:
        selected = [r for r in rows if r['case'] in
                    ([f'hall_h{helicity}_n{n}' for n in (32, 64, 128)] if kind == 'hall'
                     else [f'vacuum_n{n}' for n in (32, 64, 128)])]
        selected.sort(key=lambda r:r['cells'])
        if len(selected) != 3:
            raise RuntimeError(f'Missing refinement cases for {kind}, {helicity}')
        key = 'relative_magnetic_error' if kind == 'hall' else 'By_relative_error'
        errors = [r[key] for r in selected]
        convergence[f'{kind}_{helicity}'] = {'cells': [r['cells'] for r in selected],
            'errors': errors, 'orders': np.log2(np.array(errors[:-1])/errors[1:]).tolist()}
    controls = {}
    by_name = {r['case']:r for r in rows}
    for h in (-1, 1):
        names = [f'hall_h{h}_n128'+suffix for suffix in ('', '_half_dt', '_quarter_dt')]
        states = []
        for name in names:
            with np.load(directory/'runs'/name/'wave_samples.npz') as data:
                states.append(np.array([data['Bx_final'], data['By_final']]))
        scale = by_name[names[0]]['parameters']['amplitude']*np.sqrt(128)
        differences = [float(np.linalg.norm(states[j]-states[j+1])/scale) for j in (0, 1)]
        controls[f'hall_{h}'] = {'successive_time_state_differences': differences,
            'difference_ratio': differences[0]/differences[1],
            'energy_changes': [by_name[name]['relative_endpoint_energy_change'] for name in names]}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for h, label in [(-1, 'Ion-cyclotron branch'), (1, 'Whistler branch')]:
        curve = convergence[f'hall_{h}']
        axes[0].loglog(curve['cells'], curve['errors'], 'o-', label=label)
    curve = convergence['vacuum_None']
    axes[1].loglog(curve['cells'], curve['errors'], 'o-', label='Vacuum Maxwell wave')
    for axis, title in zip(axes, ['Hybrid PIC: fixed particle step', 'Empty domain: full Maxwell']):
        axis.set(xlabel='Cells per wavelength', ylabel='Relative magnetic-field error', title=title)
        axis.grid(True, which='both', alpha=.2)
        axis.legend(fontsize=8)
    fig.suptitle('WarpX GPU verification fixtures — not FRC formation', fontsize=12)
    figure = directory/'convergence.png'
    fig.savefig(figure, dpi=180)
    plt.close(fig)
    result = {'archive_and_physical_replay_passed': True, 'cases': rows,
              'convergence': convergence, 'time_controls': controls,
              'figure_sha256': hashlib.sha256(figure.read_bytes()).hexdigest(),
              'report_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (directory/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'audited_cases': len(rows), 'convergence': convergence}, indent=2))


if __name__ == '__main__':
    run()
