"""Conservative moving state with explicitly accounted physical current energy."""
from dataclasses import dataclass
import numpy as np
from .current_thermal import CurrentThermal
from .dual_mass import DualMass
from .finite_mass_moments import SpeciesMasses
from .mhd_state import NonphysicalState
from .toroidal_state import radius_squared


@dataclass(frozen=True)
class FiniteMassParameters:
    ion_mass: float
    model: str = 'two_fluid'
    gamma_e: float = 5/3
    gamma_i: float = 5/3
    density_scale: float = 1e20
    quadrature_order: int = 5
    spatial_order: int = 2
    linear_tolerance: float = 1e-10
    max_iterations: int = 2000

    def __post_init__(self):
        if self.model not in ('two_fluid','hall') or not np.isfinite([self.ion_mass,self.gamma_e,self.gamma_i,self.density_scale]).all():
            raise ValueError('Finite physical parameters and explicit two_fluid or hall model required')
        if self.ion_mass<=0 or min(self.gamma_e,self.gamma_i)<=1 or self.density_scale<=0:
            raise ValueError('Positive physical mass/density scale and species gammas greater than one required')
        if type(self.spatial_order) is not int or self.spatial_order not in (1,2):
            raise ValueError('Explicit spatial order 1 or 2 required')

    @property
    def masses(self):
        return SpeciesMasses(self.ion_mass) if self.model=='two_fluid' else SpeciesMasses(self.ion_mass,0.)


@dataclass
class FiniteMassState:
    conserved: np.ndarray  # rho=(mi+me)n, radial momentum, axial momentum, total energy, angular momentum
    psi: np.ndarray
    g: np.ndarray
    electron_entropy: np.ndarray

    def copy(self):return FiniteMassState(*(a.copy() for a in self.arrays()))

    def arrays(self):return self.conserved,self.psi,self.g,self.electron_entropy

    def add(self,rate,dt):return FiniteMassState(*(a+dt*b for a,b in zip(self.arrays(),rate.arrays())))


@dataclass
class FiniteMassEvaluation:
    model: CurrentThermal
    primitive: np.ndarray
    current: np.ndarray
    ion_entropy: np.ndarray
    components: dict


def evaluate(mesh,state,parameters):
    u=state.conserved;nr,nz=mesh.nr,mesh.nz
    shapes=((nr,nz,5),(nr+1,nz+1),(nr,nz),(nr+1,nz))
    if any(np.shape(a)!=s or not np.isfinite(a).all() for a,s in zip(state.arrays(),shapes)):
        raise NonphysicalState('Finite moving mass/momentum/total energy, flux and species arrays required')
    scale=max(float(np.max(abs(state.psi))),np.finfo(float).tiny)
    if np.max(abs(state.psi[0]))>1e-13*scale or np.ptp(state.psi[-1])>1e-13*scale or np.max(abs(state.psi[:,0]-state.psi[:,-1]))>1e-13*scale:
        raise NonphysicalState('Regular axis, fixed wall and periodic magnetic flux required')
    rho=u[...,0]
    if np.any(rho<=0):raise NonphysicalState('Nonpositive finite-mass density; no floor applied')
    mapping=DualMass(mesh);n=mapping.density(rho)/parameters.masses.total
    model=CurrentThermal(mesh,n,parameters.ion_mass,model=parameters.model,
        gamma_e=parameters.gamma_e,gamma_i=parameters.gamma_i,density_scale=parameters.density_scale,
        quadrature_order=parameters.quadrature_order)
    if np.any(state.electron_entropy<=0):raise NonphysicalState('Nonpositive electron entropy; no repair applied')
    packed=model.geometry.pack(state.psi,state.g);current=model.current(packed)
    samples=model.geometry.sample_electric(current);r2=radius_squared(mesh)
    vr,vz=u[...,1]/rho,u[...,2]/rho;omega=u[...,4]/(rho*r2)
    bulk=.5*rho*(vr*vr+vz*vz+r2*omega*omega)
    magnetic=model.geometry.energy_density(state.psi,state.g)
    kinetic=np.sum(model.geometry.weight*model.inertia_coefficient*np.sum(samples*samples,axis=-1),axis=(2,3))/(2*mesh.volumes)
    pe=model.electron.pressure(state.electron_entropy);electron=mapping.primal_average(pe)/(parameters.gamma_e-1)
    ion=u[...,3]-bulk-magnetic-kinetic-electron;pi=(parameters.gamma_i-1)*ion
    if not np.isfinite(pi).all() or np.any(pi<=0):
        raise NonphysicalState(f'Nonpositive finite-mass ion pressure (minimum {pi.min():.6g}); no energy repair applied')
    ui=model.ion.from_pressure(mapping.density(pi));br,bz=mesh.cell_fields(state.psi)
    q=np.stack((rho,vr,vz,pi+mapping.primal_average(pe),br,bz,omega,state.g),axis=-1)
    return FiniteMassEvaluation(model,q,current,ui,dict(bulk=bulk,magnetic=magnetic,current=kinetic,electron=electron,ion=ion))


def from_primitive(mesh,parameters,psi,g,rho,pe,pi,*,vr=0.,vz=0.,omega=0.):
    """Initialize this explicit discrete energy convention once; never reset a step."""
    shape=(mesh.nr,mesh.nz)
    rho,pi,vr,vz,omega,g=[np.broadcast_to(np.asarray(a,float),shape) for a in (rho,pi,vr,vz,omega,g)]
    mapping=DualMass(mesh);n=mapping.density(rho)/parameters.masses.total
    model=CurrentThermal(mesh,n,parameters.ion_mass,model=parameters.model,gamma_e=parameters.gamma_e,
        gamma_i=parameters.gamma_i,density_scale=parameters.density_scale,quadrature_order=parameters.quadrature_order)
    ue=model.electron.from_pressure(np.broadcast_to(pe,n.shape));ui=model.ion.from_pressure(mapping.density(pi))
    energy=model.energy_density(model.geometry.pack(psi,g),ue,ui)
    r2=radius_squared(mesh);bulk=.5*rho*(vr*vr+vz*vz+r2*omega*omega)
    total=bulk+energy['magnetic']+energy['current_kinetic']+energy['electron']+pi/(parameters.gamma_i-1)
    state=FiniteMassState(np.stack((rho,rho*vr,rho*vz,total,rho*r2*omega),axis=-1),np.asarray(psi).copy(),g.copy(),ue.copy())
    evaluate(mesh,state,parameters)
    return state
