"""Explicit nonlinear toroidal MHD verification configurations."""
from dataclasses import dataclass, asdict
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd import MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.toroidal_mhd import ToroidalMHD
from .toroidal_references import equilibrium, rotating_contact, torsional_wave, toroidal_decay


@dataclass(frozen=True)
class ToroidalCase:
    kind: str
    nr: int = 8
    nz: int = 4
    radius: float = 1.
    half_length: float = 1.
    gamma: float = 5/3
    cfl: float = .25
    order: int = 2
    end_time: float = .2
    diffusivity: float = 0.
    omega: float = .4
    current_slope: float = .3
    contrast: float = .5
    velocity: float = .8
    amplitude: float = 1e-4
    axial_mode: int = 1
    dt_max: float | None = None
    control: str = 'base'

    @property
    def name(self):
        return f'{self.kind}_nr{self.nr}_nz{self.nz}_{self.control}'

    def configuration(self):
        return {**asdict(self), 'resistivity_ohm_m': self.diffusivity*MU0,
                'scope': 'bounded compressible massless-electron MHD; periodic z, conducting slip radial wall; no vacuum/Hall coupling',
                'units': 'SI; current_slope means Btheta=sqrt(mu0)*current_slope*r; wave amplitude rad/s; decay amplitude T'}

    def reference(self, mesh, time=0.):
        if self.kind == 'equilibrium':
            return equilibrium(mesh, omega=self.omega, current_slope=self.current_slope, gamma=self.gamma)
        if self.kind == 'contact':
            return rotating_contact(mesh, time, contrast=self.contrast, omega=self.omega,
                                    current_slope=self.current_slope, velocity=self.velocity, gamma=self.gamma)
        if self.kind == 'wave':
            return torsional_wave(mesh, time, amplitude=self.amplitude, axial_mode=self.axial_mode, gamma=self.gamma)
        if self.kind == 'decay':
            return toroidal_decay(mesh, time, diffusivity=self.diffusivity, amplitude=self.amplitude, gamma=self.gamma)[0]
        raise ValueError('Unknown toroidal case kind')

    def build(self):
        mesh = MHDMesh(self.radius, self.half_length, self.nr, self.nz)
        par = MHDParameters(self.gamma, self.diffusivity*MU0, self.cfl, self.order)
        return ToroidalMHD(mesh, self.reference(mesh), par)


def protocol():
    cases = []
    for omega, current, tag in ((.7, 0., 'rotation'), (0., .5, 'pinch'), (.4, .3, 'combined')):
        for nr in (8, 16, 32):
            cases.append(ToroidalCase('equilibrium', nr=nr, omega=omega, current_slope=current, control=tag))
    for nz in (24, 48, 96, 192):
        cases.append(ToroidalCase('contact', nz=nz, omega=.7, current_slope=.4, end_time=.4))
        cases.append(ToroidalCase('wave', nz=nz, end_time=.4))
    cases += [ToroidalCase('contact', nr=16, nz=96, omega=.7, current_slope=.4, end_time=.4, control='radial16'),
              ToroidalCase('contact', nz=96, omega=.7, current_slope=.4, contrast=.8, end_time=.4, control='contrast08'),
              ToroidalCase('contact', nz=96, omega=.7, current_slope=.4, order=1, end_time=.4, control='first_order')]
    for amplitude, tag in ((1e-3, 'amplitude1e-3'), (1e-5, 'amplitude1e-5')):
        cases.append(ToroidalCase('wave', nz=96, amplitude=amplitude, end_time=.4, control=tag))
    cases += [ToroidalCase('wave', nz=192, axial_mode=3, end_time=.2, control='mode3'),
              ToroidalCase('wave', nz=96, order=1, end_time=.4, control='first_order')]
    for nr in (8, 16, 32, 64, 128):
        cases.append(ToroidalCase('decay', nr=nr, diffusivity=.1, amplitude=1e-7))
    for dt in (.002, .001, .0005, .0000625):
        cases.append(ToroidalCase('contact', nz=32, omega=.7, current_slope=.4, end_time=.05,
                                  dt_max=dt, control=f'time{dt:g}'))
    return cases
