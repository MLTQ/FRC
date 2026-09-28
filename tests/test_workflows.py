from dataclasses import replace
import json
import numpy as np
import pytest
from numpy.testing import assert_allclose
from plasma_bench.config import Experiment, Material
from plasma_bench.geometry import Grid
from plasma_bench.experiments import comparison_config
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.simulation import run
from plasma_bench.results import Result
from plasma_bench.diagnostics import characteristic_scales


def small_config():
    return Experiment(grid=Grid(.15, .3, 31, 61), duration=10e-6, dt=1e-6, output_stride=2)


def test_configuration_roundtrip_and_validation(tmp_path):
    c = small_config()
    c.save(tmp_path/"config.json")
    assert Experiment.load(tmp_path/"config.json") == c
    with pytest.raises(ValueError):
        replace(c, dt=-1)
    with pytest.raises(ValueError):
        replace(c, bias=np.nan)
    with pytest.raises(ValueError):
        Material(conductivity=-1)


def test_normalization_and_explicit_reference_scales():
    for name in ("unipolar", "rlc", "burst"):
        c = comparison_config(name)
        source = VacuumField(Grid(.02, .02, 5, 5), c.coils)
        peak_current = max(abs(c.coils[0].waveform.value(np.linspace(0, c.duration, 30001))))
        current = c.coils[0].waveform.value(1e-6)
        b_per_amp = source.at(1e-6).bz[0, 2]/current
        assert_allclose(peak_current*b_per_amp, .015, rtol=1e-6)
    scales, warnings = characteristic_scales(small_config())
    assert scales["ion_skin_depth_m"] > 0
    assert scales["magnetic_reynolds_number"] == 0
    assert scales["physical_mean_free_path_m"] is None
    assert any("Hall" in w for w in warnings)


def test_results_save_load_reproduce_and_energy_budget(tmp_path):
    c = small_config()
    first = run(c)
    second = run(c)
    first.save(tmp_path)
    loaded = Result.load(tmp_path)
    for key in first.arrays:
        assert_allclose(first.arrays[key], second.arrays[key], rtol=0, atol=0)
        assert_allclose(first.arrays[key], loaded.arrays[key], rtol=0, atol=0)
    assert loaded.metadata["configuration"] == c.to_dict()
    assert loaded.metadata["summary"]["energy_balance_relative_error"] < 1e-11
    json.dumps(loaded.metadata, allow_nan=False)
    assert first.arrays["psi"].shape[0] == len(first.metadata["topology"])


def test_zero_bias_has_undefined_reversal_and_stage0_has_no_current():
    c = replace(small_config(), bias=0, stage=0)
    result = run(c)
    assert result.metadata["summary"]["field_reversal"] is None
    assert result.metadata["summary"]["min_bz_over_bias"] is None
    assert np.all(result.arrays["j_theta"] == 0)
    json.dumps(result.metadata, allow_nan=False)


def test_reflection_symmetry_with_grid_aligned_interfaces():
    result = run(small_config(), topology=False)
    assert_allclose(result.arrays["bz"], result.arrays["bz"][:, :, ::-1], atol=1e-13)
    assert_allclose(result.arrays["j_theta"], result.arrays["j_theta"][:, :, ::-1], atol=1e-5)


def test_negative_bias_reversal_convention_and_exact_step_count():
    base = small_config()
    config = replace(base, bias=-base.bias,
                     coils=tuple(replace(coil, polarity=-coil.polarity) for coil in base.coils))
    positive, negative = run(base), run(config)
    assert_allclose(positive.arrays["min_bz_over_bias"], negative.arrays["min_bz_over_bias"])
    assert positive.metadata["summary"]["axis_field_reversal"] == negative.metadata["summary"]["axis_field_reversal"]
    assert len(Experiment().times) == 401
