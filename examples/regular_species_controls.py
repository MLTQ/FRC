"""Frozen spectra and discontinuous-state controls for species transport."""
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.linalg import expm
from plasma_bench.regular_species import RegularSpecies
from .regular_species_reference import initialize


def run_case(spec):
    _, d, u, m, field = initialize(spec)
    current = field.current(m)
    root = np.sqrt(d.volume*u[..., 0])
    size = root.size
    rng = np.random.default_rng(2937)
    discontinuous = u.copy()
    discontinuous[..., 4:6] += u[..., 0, None]*rng.uniform(-.2, .2, u[..., 4:6].shape)
    pair = RegularSpecies(m).pair(discontinuous, current)
    thermo = d.fields(discontinuous)
    budgets = []
    for k in range(2):
        rate = pair['fluid_rate'][..., k+4]
        work = d.volume*thermo['specific_internal'][..., k]*rate
        q = d.volume*thermo['entropy'][..., k]*rate
        budgets.append(dict(entropy_backward=float(abs(np.sum(d.volume*rate))/np.sum(d.volume*abs(rate))),
                            quadratic_backward=float(abs(q.sum())/abs(q).sum()),
                            pressure_work_backward=float(abs(work.sum()-current@pair['species_forcing'][k])/abs(work).sum())))
    results = {}
    arrays = dict(state=discontinuous, current=current, fluid_rate=pair['fluid_rate'],
                  species_forcing=pair['species_forcing'])
    for label, transfer, skew in (('point_unskew', 'point', False), ('point_skew', 'point', True),
                                  ('projection_skew', 'projection', True)):
        operator = RegularSpecies(m, transfer=transfer)
        def action(value):
            return m.beta/E*operator.transport(current, value/root, skew=skew)/root
        matrix = np.column_stack([action(np.eye(size)[i].reshape(root.shape)).ravel() for i in range(size)])
        eigenvalues = np.linalg.eigvals(matrix)
        spectral_radius = float(np.max(abs(eigenvalues)))
        duration = 20/spectral_radius
        initial = (root*thermo['entropy'][..., 0]).ravel()
        endpoint = expm(duration*matrix)@initial
        row = dict(max_real_over_spectral_radius=float(np.max(eigenvalues.real)/spectral_radius),
                   symmetric_part_relative=float(np.linalg.norm(matrix+matrix.T)/np.linalg.norm(matrix)),
                   constant_rate_scaled=float(np.linalg.norm(action(root))/np.linalg.norm(matrix)),
                   exponential_quadratic_drift=float((endpoint@endpoint-initial@initial)/(initial@initial)),
                   spectral_radius=spectral_radius, exponential_duration=duration)
        results[label] = row
        arrays[label+'_operator'] = matrix
        arrays[label+'_eigenvalues'] = eigenvalues
        arrays[label+'_endpoint'] = endpoint
    return dict(spec=spec, status='complete', discontinuous_budgets=budgets, spectra=results), arrays
