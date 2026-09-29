"""Independent published-unit coefficients and local closure domain checks."""
import json
import numpy as np
import pytest
from scipy.constants import atomic_mass, proton_mass, electron_mass
from plasma_bench.coulomb_transport import coulomb_coefficients,regime
from plasma_bench.config import Experiment
from plasma_bench.diagnostics import characteristic_scales


def test_collision_coefficients_match_independent_nrl_cgs_numbers():
    ne=1e20;te=5.;ti=1.;mi=39.948*atomic_mass
    a=coulomb_coefficients(ne,te,ti,mi)
    # NRL's rounded independent cgs coefficient, ne converted to cm^-3.
    assert a['electron_ion_momentum_rate_s']==pytest.approx(2.91e-6*1e14*a['coulomb_log_ei']/te**1.5,rel=.002)
    assert a['ion_ion_collision_rate_s']==pytest.approx(4.80e-8*1e14*a['coulomb_log_ii']/np.sqrt(mi/proton_mass),rel=.002)
    assert a['momentum_resistivity_ohm_m']==pytest.approx(1.03e-4*a['coulomb_log_ei']/te**1.5,rel=.002)
    assert a['coulomb_log_ei']==pytest.approx(9.29606121769283,rel=1e-12)
    b=coulomb_coefficients(ne,100.,ti,mi)
    assert b['coulomb_log_ei']==pytest.approx(12.48707453502977,rel=1e-12)


def test_regime_lengths_match_plasma_frequency_and_published_units():
    r=regime(1e20,5.,1.,39.948*atomic_mass,field=.005,current_density=1e6,length=.035,timescale=1e-5)
    assert r['electron_inertial_length_m']==pytest.approx(5.31e5/np.sqrt(1e14)/100,rel=.002)
    assert r['ion_inertial_length_m']/r['electron_inertial_length_m']==pytest.approx(np.sqrt(39.948*atomic_mass/electron_mass))
    assert r['hall_length_ratio']>4
    assert r['current_drift_m_s']==pytest.approx(62415.09074460763)
    assert not r['validity_established'] and any('Hall scale' in w for w in r['warnings'])
    json.dumps(r,allow_nan=False)


def test_invalid_classical_closure_is_unavailable_without_a_log_floor():
    r=regime(1e30,.001,.001,39.948*atomic_mass,ionization=.2)
    assert r['classical_coefficients'] is None
    assert any('no floor' in w for w in r['warnings'])
    assert any('Neutral' in w for w in r['warnings'])
    with pytest.raises(ValueError):regime(0.,1.,1.,proton_mass)
    with pytest.raises(ValueError):coulomb_coefficients(1e20,.001,100.,proton_mass)
    with pytest.raises(ValueError):coulomb_coefficients(1e20,1.,1.,proton_mass,True)


def test_existing_experiment_records_conditional_estimates_without_changing_sigma():
    config=Experiment();sigma=config.material.conductivity
    scales,warnings=characteristic_scales(config)
    assert scales['coulomb_reference']['classical_coefficients'] is not None
    assert scales['physical_collision_frequency_s'] is None
    assert config.material.conductivity==sigma
    assert any('Coulomb reference:' in w for w in warnings)
    json.dumps(scales,allow_nan=False)
