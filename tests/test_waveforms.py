import numpy as np
import pytest
from numpy.testing import assert_allclose
from plasma_bench.waveforms import DC, Gaussian, Unipolar, DampedSine, SineBurst, Sampled, Callback, to_dict, from_dict


@pytest.mark.parametrize("wave", [DC(3), Gaussian(4, .4, .2), Unipolar(5, .1, .2, .2, .3), DampedSine(2, 3, .4, .1), SineBurst(2, 3, 4, .1), Sampled((0, .3, .7, 1), (0, 2, -1, 0))])
def test_derivative_and_serialization(wave):
    t = np.array([-.123, .137, .267, .419, .781, 1.653])
    h = 1e-6
    reference = (wave.value(t+h)-wave.value(t-h))/(2*h)
    assert_allclose(wave.derivative(t), reference, rtol=1e-7, atol=2e-8)
    assert_allclose(from_dict(to_dict(wave)).value(t), wave.value(t))


def test_finite_support_and_sample_convention():
    u = Unipolar(10, 1, 1, 1, 1)
    assert_allclose(u.value([0, 1, 2, 3, 4, 5]), [0, 0, 10, 10, 0, 0], atol=1e-14)
    b = SineBurst(4, 2, 3)
    assert_allclose(b.value([-1, 0, 1.5, 2]), 0, atol=1e-25)
    s = Sampled((0, 1, 2), (0, 2, 1))
    assert_allclose(s.derivative([-1, 0, 1, 2, 3]), [0, 2, -1, 0, 0])
    assert_allclose(s.value([-1, 3]), [0, 1])


def test_invalid_parameters_and_callbacks():
    with pytest.raises(ValueError):
        Gaussian(1, 0, 0)
    with pytest.raises(ValueError):
        Sampled((0, 0), (0, 1))
    with pytest.raises(ValueError):
        to_dict(Callback(lambda t: t, lambda t: 1, "linear"))
