"""Audit completed compatible Hall studies, compare trajectories and render plots."""
from pathlib import Path
import hashlib
import numpy as np
from .compatible_hall_audit import load_all, legacy_axis
from .compatible_hall_compare import comparisons, wave_modes
from .compatible_hall_figures import figures
from .compatible_hall_records import provenance
from .verification_records import write_record


def report():
    records, archives, meshes, audit = load_all()
    legacy = legacy_axis()
    comparison = comparisons(records, archives, meshes)
    paths = figures(Path('outputs/compatible_hall_verification'), records, comparison, legacy)
    keys = ('max_total_energy_j_relative_drift', 'max_mass_kg_relative_drift',
            'max_angular_momentum_absolute_drift', 'replayed_work_relative', 'local_rate_replay_relative')
    maxima = {key: max(v.get(key, 0.) for v in audit.values()) for key in keys}
    report_sources = {f'examples/{name}.py': hashlib.sha256(Path(f'examples/{name}.py').read_bytes()).hexdigest()
                      for name in ('compatible_hall_audit', 'compatible_hall_compare', 'compatible_hall_figures', 'verify_compatible_hall')}
    # Change comparisons retain a meaningful scale for zero-mean wave quantities.
    tolerance = {}
    ref = archives['wave_tolerance_1e-12']
    from .compatible_hall_compare import state_difference
    for name in ('wave_tolerance_1e-08', 'wave_b2_n16'):
        a = archives[name]
        tolerance[name] = state_difference(meshes[name], (a['final_conserved'], a['final_psi']),
                                          (ref['final_conserved'], ref['final_psi']),
                                          (ref['initial_conserved'], ref['initial_psi']))
    result = dict(completed=True, measured_cases=len(records), spatial_cases=23, evolved_cases=40,
                  provenance=provenance(), report_sources=report_sources, audit=audit, maxima=maxima,
                  comparisons=comparison, wave_modes=wave_modes(records, archives, meshes),
                  nonlinear_tolerance_comparisons=tolerance, historical_axis=legacy,
                  figures={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  scope='Bounded conducting, periodic compressible Hall reference; cold massless electrons, adiabatic represented gas, prescribed resistance. No formation, vacuum, separate species thermal closure or finite electron inertia.')
    write_record(Path('outputs/compatible_hall_verification/report.json'), result)
    print('Audited', len(records), 'complete cases; engine', provenance()['source_sha256'])
    print('Budget/replay maxima', maxima)
    print('Tolerance comparisons', tolerance)
    print('Figures', [str(p) for p in paths])


if __name__ == '__main__':
    report()
