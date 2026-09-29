"""Frozen matrix of bounded cylindrical Hall verification cases."""
from dataclasses import dataclass, asdict
from scipy.constants import proton_mass as MP
from plasma_bench.hall_parameters import HallParameters
from plasma_bench.hall_cylinder import HallCylinderMesh, LinearHallCylinder
from .hall_wave_reference import CylindricalWave


@dataclass(frozen=True)
class HallCase:
    model: str = 'two_fluid'
    helicity: int = 1
    nr: int = 16
    nz: int = 9
    radius: float = 1.
    length: float = 2.
    number_density: float = 1e18
    ion_mass: float = MP
    bias: float = .001
    resistivity: float = 3e-5
    radial_mode: int = 1
    axial_mode: int = 1
    amplitude: float = 1e-10
    steps: int = 512
    phase_duration: float = 2.
    quadrature: int = 8
    control: str = 'base'

    @property
    def name(self):
        return f'{self.model}_h{self.helicity:+d}_nr{self.nr}_{self.control}'

    def configuration(self):
        return {**asdict(self), 'scope': 'linear incompressible, uniform plasma, fixed homogeneous cylinder traces; no exterior',
                'units': 'SI; amplitude is vector potential T m; phase_duration in radians'}

    def build(self):
        p = HallParameters(self.number_density, self.ion_mass, self.bias, self.resistivity, self.model)
        mesh = HallCylinderMesh(self.radius, self.length, self.nr, self.nz, self.quadrature)
        solver = LinearHallCylinder(mesh, p)
        reference = CylindricalWave(self.radius, self.length, self.number_density, self.ion_mass,
                                   self.bias, self.resistivity, self.helicity, self.model,
                                   self.radial_mode, self.axial_mode, self.amplitude)
        values = [v.real.copy() for v in reference.potentials(mesh.r[:, None], mesh.z[None, :])]
        for value in values:
            value[[0, -1]] = 0.
        return solver, reference, solver.from_nodal(*values)


def protocol():
    cases = []
    for model in ('mhd', 'hall', 'two_fluid'):
        radius = .01 if model == 'two_fluid' else 1.
        for h in (-1, 1):
            for nr in (8, 16, 32, 64):
                cases.append(HallCase(model=model, helicity=h, nr=nr, radius=radius, length=2*radius))
            cases.append(HallCase(model=model, helicity=h, nr=32, radius=radius, length=2*radius,
                                  resistivity=0., control='ideal'))
    # Separate temporal, axial sampling, higher radial/axial mode and integration controls.
    for h in (-1, 1):
        for steps in ((128, 256, 1024, 2048, 4096, 8192) if h == -1 else (128, 256, 1024)):
            cases.append(HallCase(helicity=h, nr=32, radius=.01, length=.02, steps=steps, control=f'time{steps}'))
        cases.append(HallCase(model='hall', helicity=h, nr=32, nz=17, control='axial17'))
        cases.append(HallCase(model='hall', helicity=h, nr=32, quadrature=12, control='quad12'))
        for nr in (16, 32, 64):
            cases.append(HallCase(model='hall', helicity=h, nr=nr, nz=17, radial_mode=2, axial_mode=3, control='higher_mode'))
    return cases
