"""Explicit consistent initial states for coupled warm Hall verification."""
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd import MHDParameters
from plasma_bench.dual_mass import DualMass
from plasma_bench.toroidal_state import from_primitive
from plasma_bench.warm_hall import WarmHallMHD
from .warm_hall_reference import WarmHallWave
from .electron_reference import ElectronReference
from .compatible_hall_references import wave_state


def initialize(case):
    wave = case['kind'] == 'wave'
    ref = (WarmHallWave if wave else ElectronReference)(**case['reference'])
    mesh = MHDMesh(ref.radius, ref.length/2, case['nr'], case['nz'])
    mapping = DualMass(mesh)
    r, z = mapping.geometry.r, mapping.geometry.z
    if wave:
        initial = wave_state(mesh, ref)
        pe = ref.electron_pressure+ref.fields(r, z)['electron_pressure']
        mi = ref.ion_mass
    else:
        cfg, k = case['initial'], 2*np.pi/ref.length
        mi = cfg['ion_mass']
        n = ref.cell_average(mesh, 'density')*(1+cfg['wall_density_contrast']*np.cos(k*mesh.z[None]))
        pe = cfg['electron_pressure']*(1+cfg['pressure_contrast']*(1-(r/ref.radius)**2)**2*np.sin(k*z)
                                       +cfg['wall_pressure_contrast']*np.cos(k*z+.31))
        rc, zc = mesh.corners
        psi = .5*cfg['bias']*rc*rc+cfg['poloidal_amplitude']*rc*rc*(1-(rc/ref.radius)**2)**3*np.cos(k*zc)
        initial = from_primitive(mesh, psi, mi*n, mapping.primal_average(pe)+cfg['ion_pressure'],
                                  toroidal_over_radius=ref.cell_average(mesh, 'g'), gamma=ref.gamma,
                                  magnetic_quadrature='compatible')
    model = WarmHallMHD(mesh, initial, pe, mi, MHDParameters(gamma=ref.gamma, cfl=case['fluid_cfl'],
                        spatial_order=case['fluid_order'], positivity_reconstruction=case['positivity_reconstruction']),
                        quadrature_order=case['quadrature'], nonlinear_tolerance=case['tolerance'],
                        hall_iterations=case['hall_iterations'], pressure_iterations=case['pressure_iterations'])
    return ref, mesh, model
