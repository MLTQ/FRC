"""Planar species fluids with independent material supports and exact vacuum."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from scipy.linalg import solve_banded


@dataclass
class SpeciesMaterialState:
    faces: np.ndarray
    velocity: np.ndarray
    thermal: np.ndarray

    def added(self, rate, scale):
        return SpeciesMaterialState(*(getattr(self, k)+scale*getattr(rate, k)
                                     for k in self.__dataclass_fields__))


class MaterialSpecies:
    """Fixed species cell numbers, independent moving faces, three momenta."""
    gamma = 5/3

    def __init__(self, numbers, masses, charges, bounds, *, left_electric=0., cfl=.15,
                 maximum_plasma_phase=.04, numerical_viscosity=(0., 0.)):
        self.numbers = np.asarray(numbers, float).copy()
        self.mass = np.asarray(masses, float); self.charge = np.asarray(charges, float)
        self.bounds = np.asarray(bounds, float); self.left_electric = float(left_electric)
        self.cfl = float(cfl); self.maximum_plasma_phase = float(maximum_plasma_phase)
        self.viscosity = np.asarray(numerical_viscosity, float)
        if self.numbers.ndim != 2 or not np.isfinite(self.numbers).all() or np.any(self.numbers <= 0):
            raise ValueError('Positive finite cell particle amounts; omit empty species')
        self.species, self.cells = self.numbers.shape
        if (self.mass.shape != (self.species,) or self.charge.shape != self.mass.shape or
                not np.isfinite(self.mass+self.charge).all() or np.any(self.mass <= 0)):
            raise ValueError('Finite species charges and physical positive masses required')
        if (self.bounds.shape != (2,) or not np.isfinite(self.bounds).all() or self.bounds[0] >= self.bounds[1] or
                not np.isfinite([left_electric,cfl,maximum_plasma_phase]).all() or cfl <= 0 or maximum_plasma_phase <= 0 or
                self.viscosity.shape != (2,) or not np.isfinite(self.viscosity).all() or np.any(self.viscosity < 0)):
            raise ValueError('Ordered remote bounds and positive finite time controls required')
        cell_mass = self.mass[:, None]*self.numbers
        self.matrix = np.zeros((self.species, 3, self.cells+1))
        self.matrix[:, 1, :-1] += cell_mass/3; self.matrix[:, 1, 1:] += cell_mass/3
        self.matrix[:, 0, 1:] = cell_mass/6; self.matrix[:, 2, :-1] = cell_mass/6

    def validate(self, state):
        if (state.faces.shape != (self.species,self.cells+1) or state.velocity.shape != (self.species,self.cells+1,3)
                or state.thermal.shape != self.numbers.shape or not all(np.isfinite(getattr(state,k)).all() for k in state.__dataclass_fields__)):
            raise ValueError('Complete finite material species state required')
        width = np.diff(state.faces, axis=1)
        if np.any(width <= 0) or np.any(state.faces[:,0] <= self.bounds[0]) or np.any(state.faces[:,-1] >= self.bounds[1]):
            raise ValueError('Species cells crossed or reached a remote field boundary; no remapping')
        if np.any(state.thermal < 0):
            raise ValueError('Negative thermal energy; no reset or floor')
        return width

    def state(self, faces, temperature_ev, velocity=None):
        faces = np.asarray(faces, float).copy()
        velocity = np.zeros((self.species,self.cells+1,3)) if velocity is None else np.asarray(velocity,float).copy()
        thermal = self.numbers*EV*np.broadcast_to(temperature_ev,self.numbers.shape)/(self.gamma-1)
        result = SpeciesMaterialState(faces,velocity,thermal.copy()); self.validate(result); return result

    def field(self, state):
        width = self.validate(state)
        edges = np.unique(np.r_[self.bounds,state.faces.ravel()]); h = np.diff(edges); centre = .5*(edges[:-1]+edges[1:])
        charge = np.zeros_like(centre); active = []
        for s in range(self.species):
            cell = np.searchsorted(state.faces[s],centre,side='right')-1
            inside = (cell >= 0)&(cell < self.cells); active.append((cell,inside))
            charge[inside] += EV*self.charge[s]*self.numbers[s,cell[inside]]/width[s,cell[inside]]
        electric = np.r_[self.left_electric,self.left_electric+np.cumsum(charge*h)/EPS0]
        energy = EPS0/6*np.sum(h*(electric[:-1]**2+electric[:-1]*electric[1:]+electric[1:]**2))
        return edges,electric,active,float(energy)

    def electric_force(self, state, geometry=None):
        edges,electric,active,_ = self.field(state) if geometry is None else geometry
        h = np.diff(edges); node = np.array([.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)])
        points = edges[:-1,None]+h[:,None]*node
        field = electric[:-1,None]+np.diff(electric)[:,None]*node
        force = np.zeros((self.species,self.cells+1))
        for s,(cell,inside) in enumerate(active):
            j=cell[inside]; width=np.diff(state.faces[s])[j]
            right=(points[inside]-state.faces[s,j,None])/width[:,None]
            coefficient=EV*self.charge[s]*self.numbers[s,j]/width*h[inside]/2
            np.add.at(force[s],j,coefficient*np.sum(field[inside]*(1-right),axis=1))
            np.add.at(force[s],j+1,coefficient*np.sum(field[inside]*right,axis=1))
        return force

    def pressure(self, state):
        width = self.validate(state); pressure=(self.gamma-1)*state.thermal/width
        sound=np.sqrt(self.gamma*pressure*width/(self.mass[:,None]*self.numbers))
        compression=np.maximum(-np.diff(state.velocity[...,2],axis=1),0.)
        viscosity=self.mass[:,None]*self.numbers/width*(self.viscosity[0]*sound*compression+self.viscosity[1]*compression**2)
        return pressure,viscosity,sound

    def derivative(self, state):
        p,q,_=self.pressure(state); stress=p+q; force=np.zeros_like(state.velocity)
        force[:,:-1,2]-=stress;force[:,1:,2]+=stress
        force[...,2]+=self.electric_force(state)
        acceleration=np.stack([solve_banded((1,1),m,f,check_finite=False) for m,f in zip(self.matrix,force)])
        return SpeciesMaterialState(state.velocity[...,2],acceleration,-stress*np.diff(state.velocity[...,2],axis=1))

    def timestep(self, state):
        width=self.validate(state);_,_,sound=self.pressure(state); gradient=abs(np.diff(state.velocity[...,2],axis=1))
        speed=sound+gradient
        acoustic=np.min(np.divide(self.cfl*width,speed,out=np.full_like(width,np.inf),where=speed>0))
        # Sum of species maximum densities bounds the local overlap plasma rate.
        omega2=EV**2/EPS0*np.sum(self.charge**2/self.mass*np.max(self.numbers/width,axis=1))
        return min(acoustic,self.maximum_plasma_phase/np.sqrt(omega2)) if omega2 else acoustic

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite material step required')
        k1=self.derivative(state);k2=self.derivative(state.added(k1,dt/2))
        k3=self.derivative(state.added(k2,dt/2));k4=self.derivative(state.added(k3,dt))
        rate=SpeciesMaterialState(*(sum(w*getattr(k,name) for w,k in zip((1,2,2,1),(k1,k2,k3,k4)))/6 for name in state.__dataclass_fields__))
        result=state.added(rate,dt);self.validate(result);return result

    def sample(self,state,points):
        points=np.asarray(points,float); width=self.validate(state)
        if points.ndim!=1 or not np.isfinite(points).all():raise ValueError('Finite one-dimensional sampling coordinates required')
        number=np.zeros((self.species,)+points.shape); velocity=np.zeros(number.shape+(3,)); temperature=np.zeros_like(number)
        for s in range(self.species):
            j=np.searchsorted(state.faces[s],points,side='right')-1;inside=(j>=0)&(j<self.cells);k=j[inside]
            a=(points[inside]-state.faces[s,k])/width[s,k]
            number[s,inside]=self.numbers[s,k]/width[s,k]
            velocity[s,inside]=(1-a[:,None])*state.velocity[s,k]+a[:,None]*state.velocity[s,k+1]
            temperature[s,inside]=(self.gamma-1)*state.thermal[s,k]/(self.numbers[s,k]*EV)
        edges,electric,_,_=self.field(state)
        current=EV*np.sum(self.charge[:,None,None]*number[...,None]*velocity,axis=0)
        return dict(number=number,velocity=velocity,temperature_ev=temperature,current=current,electric=np.interp(points,edges,electric))

    def totals(self,state):
        self.validate(state);left=state.velocity[:,:-1];right=state.velocity[:,1:];mass=self.mass[:,None]*self.numbers
        kinetic=np.sum(mass/6*np.sum(left*left+left*right+right*right,axis=-1))
        momentum=np.sum(mass[...,None]*(left+right)/2,axis=(0,1))
        edges,electric,active,energy=self.field(state)
        empty=~np.any(np.stack([inside for _,inside in active]),axis=0)
        exterior=EPS0/6*np.sum(np.diff(edges)[empty]*(electric[:-1][empty]**2+electric[:-1][empty]*electric[1:][empty]+electric[1:][empty]**2))
        return dict(numbers=self.numbers.sum(axis=1),charge_c_m2=float(EV*np.sum(self.charge*self.numbers.sum(axis=1))),
                    momentum_kg_per_m_s=momentum,kinetic_j_m2=float(kinetic),thermal_j_m2=float(state.thermal.sum()),
                    field_j_m2=energy,empty_field_j_m2=float(exterior),energy_j_m2=float(kinetic+state.thermal.sum()+energy))
