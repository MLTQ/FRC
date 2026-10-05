"""Independent full Ohm-field curl and divergence of the continuum reference."""
import numpy as np
from scipy.special import j0, j1
from examples.continuum_contact_mode import ContinuumContactMode


def test_full_linear_electric_field_satisfies_continuum_faraday():
    mode = ContinuumContactMode(); radius = mode.ref.radius/3
    z = np.arange(mode.points)*mode.ref.length/mode.points-mode.ref.length/2
    fields = mode.fields(radius, z)
    electric = (fields['electric']/np.array([j1(mode.q*radius), j1(mode.q*radius), j0(mode.q*radius)])).T.ravel()
    magnetic = mode.physical_nodes[6:].ravel()
    assert np.linalg.norm(-mode.curl@electric-mode.growth*magnetic)/np.linalg.norm(mode.growth*magnetic) < 5e-9
    divergence = mode.q*mode.physical_nodes[6]+mode.derivative@mode.physical_nodes[8]
    scale = np.linalg.norm(mode.q*mode.physical_nodes[6])+np.linalg.norm(mode.derivative@mode.physical_nodes[8])
    assert np.linalg.norm(divergence)/scale < 1e-8
    wall = mode.fields(mode.ref.radius, z)
    for k in (1, 6, 7):
        assert np.linalg.norm(wall['physical'][..., k])/np.linalg.norm(fields['physical'][..., k]) < 1e-14
