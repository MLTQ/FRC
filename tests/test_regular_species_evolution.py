"""Independent full-transit shear and temporal convergence controls."""
import numpy as np
from examples.regular_species_evolution import run_case, audit


def test_full_transit_shear_refines_and_preserves_actual_entropy():
    coarse_spec = dict(elements=4)
    fine_spec = dict(elements=8)
    coarse, ca = run_case(coarse_spec)
    fine, fa = run_case(fine_spec)
    c = audit(coarse_spec, ca, 12)
    f = audit(fine_spec, fa, 12)
    assert c['signal_error']['relative_L2']/f['signal_error']['relative_L2'] > 10
    assert f['signal_error']['relative_L2'] < .003
    assert abs(fine['grid_entropy_drift_per_particle']) < 2e-13
    assert abs(fine['grid_quadratic_relative_drift']) < 2e-11
    assert f['prescribed_current_error'] < 2e-12


def test_time_refinement_changes_converge_at_rk4_order():
    endpoints = [run_case(dict(elements=4, steps=n))[1]['final_entropy'] for n in (48, 96, 192)]
    ratio = np.linalg.norm(endpoints[0]-endpoints[1])/np.linalg.norm(endpoints[1]-endpoints[2])
    assert 14 < ratio < 18
