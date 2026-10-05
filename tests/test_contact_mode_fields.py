"""Independent cylindrical polynomial control of the mode's strong curl."""
import numpy as np
from examples.contact_mode_fields import magnetic_curl
from examples.current_transport_sampling import TransportSampler
from examples.regular_plasma_reference import initialize


def test_strong_curl_matches_original_cylindrical_polynomial():
    _, d, _, system, _ = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    sample = TransportSampler(d, system.geometry.mesh, 16); mesh = system.geometry.mesh
    x = .5*mesh.rf**2; edge = x[-1]
    psi = np.broadcast_to((x*(edge-x))[:, None], (mesh.nr+1, mesh.nz))
    g = np.full((mesh.nr, mesh.nz), .23)
    value = np.r_[psi.ravel(), g.ravel()]
    actual = magnetic_curl(sample, value)
    target = np.zeros_like(actual); target[..., 1] = 2*sample.R; target[..., 2] = .46
    assert np.max(abs(actual-target)) < 2e-12
