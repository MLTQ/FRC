"""Audited electron-pressure convergence and field/thermal-exchange figures."""
from pathlib import Path
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter, NullFormatter
from .electron_audit import load, comparisons, physical_drift
from .electron_records import provenance
from .verification_records import write_record


def precision_controls(records):
    """Check that the intentionally strict trial does not recover with more iterations."""
    from .electron_measurements import initialize
    from plasma_bench.electron_midpoint import step, ElectronConvergenceError
    case = records['battery_roundoff_control']['configuration']
    _, _, model, g, u = initialize(case)
    results = []
    for budget in (160, 320):
        oldg, oldu = g.copy(), u.copy()
        try:
            step(model, g, u, case['end_time']/case['steps'], tolerance=case['tolerance'], max_iterations=budget)
        except ElectronConvergenceError as error:
            results.append(dict(max_iterations=budget, failure=str(error), accepted_time_s=0.))
        else:
            raise AssertionError('Strict precision control unexpectedly converged; inspect and update interpretation')
        np.testing.assert_array_equal(g, oldg); np.testing.assert_array_equal(u, oldu)
    return results


def figures(directory, records, arrays, comparison, physical):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    paths = []
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.1), constrained_layout=True)
    nr = [8, 16, 32, 64, 128]
    for key, label in [('battery_relative_L2', 'Pressure-driven ġ'), ('thermal_rate_relative_L2', 'Electron thermal drift'),
                       ('axis_battery_relative_L2', 'First-annulus ġ')]:
        axes[0].loglog(nr, [records[f'spatial_n{n}']['metrics'][key] for n in nr], 'o-', label=label)
    axes[0].set(title='Independent analytic spatial rates', xlabel='Radial cells', ylabel='Relative L2 error')
    nz = [16, 32, 64, 128]
    axes[1].loglog(nz, [records[f'drift_nz{n}']['metrics']['pressure_perturbation_relative_L2'] for n in nz], 'o-', label='One electron transit; nodal pressure')
    axes[1].loglog(nz, [physical[f'drift_nz{n}']['relative_L2'] for n in nz], 's--', label='Continuous reconstructed pressure')
    axes[1].set(title='Actual heat translation, fixed bulk', xlabel='Axial cells', ylabel='Pressure error / perturbation')
    for family, counts in [('drift', (64, 128, 256)), ('nonlinear', (10, 20, 40))]:
        axes[2].loglog(counts, [comparison['time'][f'{family}_time_{n}']['pressure_relative_to_reference_change'] for n in counts], 'o-', label=family+' pressure')
    axes[2].set(title='Fixed-grid time convergence', xlabel='Accepted steps', ylabel='Error / fine-reference change')
    for ax, ticks in zip(axes, [nr, nz, [10, 20, 40, 80, 160, 320]]):
        ax.grid(True, which='both', alpha=.2); ax.legend(fontsize=8)
        ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    fig.suptitle('Electron pressure and entropy component • independent errors, not a full fluid model')
    p = directory/'convergence.png'; fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)

    name = 'battery_n64'; a, record = arrays[name], records[name]
    rf, zf = a['radial_faces'], a['axial_faces']
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.1), constrained_layout=True)
    bt = .5*(rf[1:]+rf[:-1])[:, None]*a['final_g']*1e3
    pc = axes[0].pcolormesh(zf*100, rf*100, bt, cmap='RdBu_r', vmin=-abs(bt).max(), vmax=abs(bt).max())
    fig.colorbar(pc, ax=axes[0], label='Bθ at cell centers [mT]')
    axes[0].set(title='Generated toroidal field at 1 μs', xlabel='z [cm]', ylabel='r [cm]')
    dp = a['final_pressure']-a['initial_pressure']
    dual_r = np.sqrt(np.r_[rf[0]**2, .5*(rf[:-1]**2+rf[1:]**2), rf[-1]**2])
    dual_z = np.r_[zf[0], .5*(zf[:-1]+zf[1:]), zf[-1]]
    pc = axes[1].pcolormesh(dual_z*100, dual_r*100, np.concatenate((dp, dp[:, :1]), axis=1), cmap='RdBu_r', vmin=-abs(dp).max(), vmax=abs(dp).max())
    axes[1].set_xlim(zf[0]*100,zf[-1]*100); axes[1].set_ylim(0,rf[-1]*100)
    fig.colorbar(pc, ax=axes[1], label='Electron pressure change [Pa]')
    axes[1].set(title='Coupled electron thermal response', xlabel='z [cm]', ylabel='r [cm]')
    history = record['history_columns']; time = np.asarray(history['time_s'])*1e6
    axes[2].plot(time, np.asarray(history['magnetic_energy_j'])*1e6, label='Magnetic energy gain')
    axes[2].plot(time, (history['electron_internal_energy_j'][0]-np.asarray(history['electron_internal_energy_j']))*1e6, '--', label='Electron thermal energy loss')
    axes[2].set(title='Resolved pressure–field exchange', xlabel='Time [μs]', ylabel='Energy [μJ]')
    axes[2].grid(True, alpha=.2); axes[2].legend(fontsize=8)
    fig.suptitle('Pressure battery from Bθ=0 • fixed positive density • no ion/fluid or Hall evolution')
    p = directory/'battery.png'; fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)
    return paths


def run():
    directory = Path('outputs/electron_pressure_verification')
    records, arrays, audits = load(directory)
    comparison = comparisons(records, arrays)
    physical = {}
    for name, record in records.items():
        if record['configuration']['kind'] == 'drift':
            first, check = physical_drift(record, arrays[name], 8), physical_drift(record, arrays[name], 12)
            physical[name] = dict(relative_L2=first, quadrature_change=abs(first-check))
    paths = figures(directory, records, arrays, comparison, physical)
    keys = ('charge_backward_error', 'work_backward_error', 'logarithmic_entropy_rate_backward_error',
            'max_energy_relative_drift', 'max_logarithmic_entropy_change_per_particle', 'max_proxy_relative_drift')
    maxima = {k: max(v.get(k, 0.) for v in audits.values()) for k in keys}
    sources = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path('examples/electron_audit.py'), Path(__file__))}
    write_record(directory/'report.json', dict(completed=True, provenance=provenance(), report_sources=sources,
                  cases=len(records), expected_failures=sum(r['configuration']['expected_failure'] for r in records.values()),
                  audits=audits, maxima=maxima, comparisons=comparison, physical_drift=physical,
                  precision_controls=precision_controls(records),
                  figures={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
    print('Audited cases', len(records), 'maxima', maxima, 'comparisons', comparison)


if __name__ == '__main__':
    run()
