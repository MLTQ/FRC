"""Fully recorded fixed-material current-response verification inputs."""
from dataclasses import dataclass,asdict
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.current_inertia import MaterialCurrent
from .inertial_sphere_reference import InertialSphere
from .inertial_screening_reference import InertialScreening
from .profiled_current_reference import ProfiledSphere


@dataclass(frozen=True)
class CurrentCase:
    name: str
    kind: str
    radial_cells: int
    electron_length: float=.2
    dt: float=.002
    end_time: float=.1
    matching_radius: float=2.
    angular_cells: int|None=None
    modes: int=8
    quadrature: int=5
    grading: float=1.

    def reference(self):
        if self.kind=='decay':return InertialSphere(electron_length=self.electron_length)
        if self.kind=='screening':return InertialScreening(electron_length=self.electron_length)
        if self.kind=='profiled':return ProfiledSphere(electron_length=self.electron_length)
        raise ValueError('Unknown current verification case')

    def build(self):
        d=MovingDomain(self.radial_cells,self.angular_cells or 2*self.radial_cells,
                       matching_radius=self.matching_radius,vacuum_layers=round(self.radial_cells*(self.matching_radius-1)),
                       modes=self.modes,quadrature=self.quadrature,radial_grading=self.grading)
        ref=self.reference()
        density=ref.density if self.kind=='profiled' else ref.electron_density
        eta=ref.resistivity if self.kind=='profiled' else 0. if self.kind=='screening' else .1*MU0
        operator=MaterialCurrent(d,d.fluid_reference,electron_density=density,resistivity=eta)
        initial=np.zeros(d.fluid_node_count) if self.kind=='screening' else ref.fields(*d.fluid_reference.T)['a']
        return operator,ref,initial

    def configuration(self):
        ref=self.reference()
        return {**asdict(self),'angular_cells':self.angular_cells or 2*self.radial_cells,
                'reference_parameters':asdict(ref),'radius_m':1.,'physical_electron_mass':True,
                'massless_comparison':self.electron_length==0,'initialization':'zero field/current' if self.kind=='screening' else 'nodal interpolation of continuum a',
                'resistivity_profile':'lambda*mu0*(k^-2+de(s)^2)' if self.kind=='profiled' else 'zero' if self.kind=='screening' else '0.1*mu0 ohm m',
                'density_profile':'me/[mu0 e^2 d0^2 (1+2 s^2/R^2)]' if self.kind=='profiled' else 'me/(mu0 e^2 de^2), or explicit massless comparison',
                'bias':'linear 0 to 1e-8 T over end_time' if self.kind=='screening' else 'zero',
                'geometry':'fixed, radius-1 sphere approximated by actual material mesh',
                'exterior':'elliptic true vacuum to infinity; no exterior carriers',
                'interpretation':'Linear current component; no moving thermal edge, Hall feedback, transport evolution or formation claim'}


def protocol():
    cases=[CurrentCase(f'decay_d{de:g}_{n}','decay',n,de) for de in (0.,.05,.2,1.) for n in (8,16,32)]
    cases += [CurrentCase('decay_d0.2_64','decay',64,.2)]
    cases += [CurrentCase(f'screen_d{de:g}_{n}','screening',n,de,dt=.01,end_time=.01) for de in (.05,.2,1.,4.) for n in (8,16,32)]
    cases += [CurrentCase('screen_d0.05_64','screening',64,.05,dt=.01,end_time=.01)]
    cases += [CurrentCase(f'profiled_{n}','profiled',n,.3) for n in (8,16,32)]
    cases += [CurrentCase(f'decay_d0.2_16_{tag}','decay',16,.2,**kwargs) for tag,kwargs in
              [('time',{'dt':.001}),('matching',{'matching_radius':3.}),('angular',{'angular_cells':64}),
               ('quadrature',{'quadrature':7}),('grading',{'grading':2.}),('modes',{'modes':12})]]
    return cases
