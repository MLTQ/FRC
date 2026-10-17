"""Complete NIST-table and Maxwellian-integral comparison for three molecules."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.constants import elementary_charge, electron_mass
from plasma_bench.molecular_beb import MolecularBEB, BOHR_RADIUS_M, RYDBERG_EV
from .implicit_plasma_records import provenance, freeze


def adaptive_integral(model, temperature):
    result = 0.
    for j, (binding, _, _) in enumerate(model.orbitals):
        upper = (model.maximum_energy_ev-binding)/temperature
        # Adaptive Gauss-Kronrod evaluates the same bounded integral using
        # independent nodes/error control. Scaling makes its tolerance explicit.
        def integrand(x):
            energy = binding+temperature*x
            return float(model.orbital_cross_sections(energy)[j])/1e-20*energy*np.exp(-x)
        value = quad(integrand, 0., upper, epsabs=1e-10, epsrel=1e-11, limit=300)[0]*1e-20
        result += np.sqrt(8*elementary_charge/(np.pi*electron_mass*temperature))*np.exp(-binding/temperature)*value
    return result


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh molecular comparison output required')
    out.mkdir(parents=True)
    identity = provenance(['examples.molecular_beb'])
    freeze(out, identity)
    root = Path('data/molecular/nist_beb')
    manifest = json.loads((root/'manifest.json').read_text())
    report = dict(terminal=False, source_sha256=identity['source_sha256'], gases=[],
        input_manifest=manifest, table_rounding_half_unit_angstrom2=.0005,
        physical_constants=dict(bohr_radius_m=BOHR_RADIUS_M, rydberg_energy_ev=RYDBERG_EV,
                                source='NIST CODATA 2022'),
        scope='Numerical reproduction of NIST BEB theoretical total-ionization tables and bounded '
              'Maxwellian integration. No molecular products, energy partition, plasma evolution or lab validation.')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    temperatures = np.array([1., 2., 5., 10., 20., 50., 100.])
    for gas in ('H2', 'N2', 'O2'):
        model = MolecularBEB.from_nist(root, gas)
        path = root/f'{gas}_cross_sections.txt'
        if hashlib.sha256(path.read_bytes()).hexdigest() != manifest['files'][path.name]['sha256']:
            raise ValueError('Reference cross sections changed since acquisition')
        table = np.loadtxt(path, skiprows=3)
        predicted = model.cross_section(table[:, 0])/1e-20
        error = predicted-table[:, 1]
        integrals = {order: [model.maxwellian_integral(t, order) for t in temperatures] for order in (8, 16, 32)}
        adaptive = np.array([adaptive_integral(model, t) for t in temperatures])
        rates = {order: np.array([r['rate_within_energy_range_m3_s'] for r in rows])
                 for order, rows in integrals.items()}
        row = dict(gas=gas, table_points=len(table), maximum_table_difference_angstrom2=float(abs(error).max()),
            points_outside_published_rounding=int(np.sum(abs(error) > .0005)),
            maximum_relative_quadrature_difference={order:float(np.max(abs(value-adaptive)/adaptive))
                for order,value in rates.items()}, integrals=integrals, adaptive_rates_m3_s=adaptive.tolist())
        report['gases'].append(row)
        np.savez_compressed(out/f'{gas}.npz', energy_ev=table[:, 0], nist_angstrom2=table[:, 1],
            predicted_angstrom2=predicted, temperature_ev=temperatures, adaptive_m3_s=adaptive,
            **{f'gauss_order_{order}_m3_s':value for order,value in rates.items()})
        energy = np.geomspace(min(model.orbitals[:, 0]), model.maximum_energy_ev, 500)
        line, = axes[0].semilogx(energy, model.cross_section(energy)/1e-20, label=gas)
        axes[0].scatter(table[:, 0], table[:, 1], s=8, marker='x', color=line.get_color())
        axes[1].loglog(temperatures, rates[32], 'o-', label=gas)
        print({k:v for k,v in row.items() if k not in ('integrals','adaptive_rates_m3_s')}, flush=True)
    axes[0].set(xlabel='Incident electron energy [eV]', ylabel='Total ionization cross section [Å²]',
                title='Orbital formula and independent NIST tables (×)')
    axes[1].set(xlabel='Electron temperature [eV]', ylabel='Rate integral, E ≤ 5000 eV [m³/s]',
                title='Stationary targets; Maxwellian electrons')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend()
    fig.savefig(out/'molecular_ionization.png', dpi=160)
    fig.savefig(out/'molecular_ionization.pdf')
    plt.close(fig)
    report['terminal'] = True
    # Report the actual table discrepancies, including any excess over its
    # published rounding, instead of fitting parameters to force agreement.
    report['success'] = all(r['points_outside_published_rounding'] == 0 and
        r['maximum_relative_quadrature_difference'][32] < 1e-9 for r in report['gases'])
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    if identity != provenance(['examples.molecular_beb']):
        raise RuntimeError('Source changed during molecular comparison')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    result = run(parser.parse_args().destination)
    raise SystemExit(0 if result['success'] else 1)
