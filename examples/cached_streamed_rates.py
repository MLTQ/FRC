"""Reuse unchanged axial curl bases inside complete streamed spatial audits."""
import numpy as np
from . import streamed_radial_rates as driver
from . import streamed_interface_audit as interface
from .current_transport_sampling import bases
from .implicit_plasma_records import provenance


def interfaces(*args, **options):
    original = interface.magnetic_curl
    axial_second = None
    def magnetic_curl(sample, magnetic):
        nonlocal axial_second
        mesh = sample.mesh; count = (mesh.nr+1)*mesh.nz
        psi = magnetic[:count].reshape(mesh.nr+1, mesh.nz)
        g = magnetic[count:].reshape(mesh.nr, mesh.nz)
        dxx = bases(.5*mesh.rf**2, .5*sample.r**2, 5)[2]
        if axial_second is None:
            axial_second = bases(mesh.zf, sample.z, 5, True)[2]
        return np.stack((-sample.R*(sample.ex@g@sample.dez.T),
            -(sample.nx@psi@axial_second.T)/sample.R-sample.R*(dxx@psi@sample.nz.T),
            (2*sample.ex+sample.R**2*sample.dex)@g@sample.ez.T), axis=-1)
    try:
        interface.magnetic_curl = magnetic_curl
        return interface.measure(*args, **options)
    finally:
        interface.magnetic_curl = original


def run(destination, specs, **options):
    original = driver.interfaces, driver.provenance
    try:
        driver.interfaces = interfaces
        driver.provenance = lambda seeds: provenance([*seeds, 'examples.cached_streamed_rates',
            'examples.streamed_interface_audit', 'examples.streamed_radial_rates'])
        return driver.run(destination, specs, **options)
    finally:
        driver.interfaces, driver.provenance = original
