"""Compatible material transport with different cell counts for each species."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from scipy.linalg import solve_banded
from .species_material import SpeciesMaterialState


@dataclass
class MaterialSetState:
    blocks: tuple

    @property
    def faces(self):
        return np.concatenate([s.faces.ravel() for s in self.blocks])

    def added(self,rate,scale):
        return MaterialSetState(tuple(s.added(r,scale) for s,r in zip(self.blocks,rate.blocks)))


class MaterialSpeciesSet:
    """One occupied material block per connected species support; no padding."""
    def __init__(self,models,species_indices,*,reference_velocity=(0.,0.,0.)):
        self.models=tuple(models);self.species_indices=np.asarray(species_indices,int)
        self.reference_velocity=np.asarray(reference_velocity,float).copy()
        if self.reference_velocity.shape!=(3,) or not np.isfinite(self.reference_velocity).all():raise ValueError('Finite fixed material velocity frame required')
        if not self.models or len(self.models)!=len(self.species_indices) or np.any(self.species_indices<0):
            raise ValueError('Nonempty material blocks with explicit atomic species IDs required')
        if any(m.species!=1 for m in self.models):raise ValueError('Every support block contains one species')
        self.bounds=self.models[0].bounds;self.left_electric=self.models[0].left_electric
        if any(not np.array_equal(m.bounds,self.bounds) or m.left_electric!=self.left_electric for m in self.models):
            raise ValueError('All species share physical remote electric boundaries')
        self.mass=np.array([m.mass[0] for m in self.models]);self.charge=np.array([m.charge[0] for m in self.models]);self.species=len(self.models)

    def validate(self,state):
        if len(state.blocks)!=self.species:raise ValueError('A state for every occupied support block required')
        widths=[m.validate(s) for m,s in zip(self.models,state.blocks)]
        for index in set(self.species_indices):
            supports=sorted((s.faces[0,0],s.faces[0,-1]) for k,s in zip(self.species_indices,state.blocks) if k==index)
            if any(a[1]>b[0] for a,b in zip(supports[:-1],supports[1:])):
                raise ValueError('Disconnected pieces of one species overlap; no implicit multistream closure')
        return widths

    def field(self,state):
        widths=self.validate(state);edges=np.unique(np.r_[self.bounds,state.faces]);h=np.diff(edges);x=(edges[:-1]+edges[1:])/2
        charge=np.zeros_like(x);active=[]
        for m,s,width in zip(self.models,state.blocks,widths):
            j=np.searchsorted(s.faces[0],x,side='right')-1;inside=(j>=0)&(j<m.cells);active.append((j,inside))
            charge[inside]+=EV*m.charge[0]*m.numbers[0,j[inside]]/width[0,j[inside]]
        electric=np.r_[self.left_electric,self.left_electric+np.cumsum(charge*h)/EPS0]
        energy=EPS0/6*np.sum(h*(electric[:-1]**2+electric[:-1]*electric[1:]+electric[1:]**2))
        return edges,electric,active,float(energy)

    def derivative(self,state):
        edges,electric,active,energy=self.field(state);rates=[]
        for m,s,occupied in zip(self.models,state.blocks,active):
            p,q,_=m.pressure(s);stress=p+q;force=np.zeros_like(s.velocity)
            force[:,:-1,2]-=stress;force[:,1:,2]+=stress
            force[...,2]+=m.electric_force(s,(edges,electric,[occupied],energy))
            acceleration=solve_banded((1,1),m.matrix[0],force[0],check_finite=False)[None,...]
            rates.append(SpeciesMaterialState(s.velocity[...,2]+self.reference_velocity[2],acceleration,-stress*np.diff(s.velocity[...,2],axis=1)))
        return MaterialSetState(tuple(rates))

    def timestep(self,state):
        widths=self.validate(state);acoustic=np.inf;omega2=0.;phase=np.inf
        for m,s,width in zip(self.models,state.blocks,widths):
            _,_,sound=m.pressure(s);speed=sound+abs(np.diff(s.velocity[...,2],axis=1))
            acoustic=min(acoustic,float(np.min(np.divide(m.cfl*width,speed,out=np.full_like(width,np.inf),where=speed>0))))
            omega2+=EV**2/EPS0*m.charge[0]**2/m.mass[0]*np.max(m.numbers/width);phase=min(phase,m.maximum_plasma_phase)
        return min(acoustic,phase/np.sqrt(omega2)) if omega2 else acoustic

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite material step required')
        k1=self.derivative(state);k2=self.derivative(state.added(k1,dt/2));k3=self.derivative(state.added(k2,dt/2));k4=self.derivative(state.added(k3,dt))
        blocks=[]
        for values in zip(k1.blocks,k2.blocks,k3.blocks,k4.blocks):
            blocks.append(SpeciesMaterialState(*(sum(w*getattr(s,name) for w,s in zip((1,2,2,1),values))/6 for name in SpeciesMaterialState.__dataclass_fields__)))
        result=state.added(MaterialSetState(tuple(blocks)),dt);self.validate(result);return result

    def sample(self,state,points):
        self.validate(state);values=[m.sample(s,points) for m,s in zip(self.models,state.blocks)]
        result={key:np.concatenate([v[key] for v in values],axis=0) for key in ('number','velocity','temperature_ev')}
        result['velocity']=np.where(result['number'][...,None]>0,result['velocity']+self.reference_velocity,0.)
        result['current']=EV*np.sum(self.charge[:,None,None]*result['number'][...,None]*result['velocity'],axis=0)
        edges,electric,_,_=self.field(state);result['electric']=np.interp(points,edges,electric);return result

    def totals(self,state):
        self.validate(state);number=[];momentum=np.zeros(3);kinetic=0.;thermal=0.
        for m,s in zip(self.models,state.blocks):
            left,right=s.velocity[0,:-1]+self.reference_velocity,s.velocity[0,1:]+self.reference_velocity;mass=m.mass[0]*m.numbers[0]
            number.append(m.numbers.sum());momentum+=np.sum(mass[:,None]*(left+right)/2,axis=0)
            kinetic+=np.sum(mass/6*np.sum(left*left+left*right+right*right,axis=1));thermal+=s.thermal.sum()
        edges,electric,active,field=self.field(state);empty=~np.any(np.stack([inside for _,inside in active]),axis=0)
        exterior=EPS0/6*np.sum(np.diff(edges)[empty]*(electric[:-1][empty]**2+electric[:-1][empty]*electric[1:][empty]+electric[1:][empty]**2))
        return dict(numbers_m2=np.asarray(number),momentum_kg_per_m_s=momentum,kinetic_j_m2=float(kinetic),thermal_j_m2=float(thermal),
                    electric_energy_j_m2=field,empty_field_j_m2=float(exterior),energy_j_m2=float(kinetic+thermal+field))
