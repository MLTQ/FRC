"""Local atomic moment evolution on the true overlap of material species."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0


@dataclass
class ReactionOverlapState:
    edges: np.ndarray
    points: np.ndarray
    weights: np.ndarray
    moments: np.ndarray
    left_electric: float


class MaterialReactionOverlap:
    """A frozen-geometry reaction substep, without lossy grid redeposition."""
    def __init__(self, source, reference_density, *, quadrature=3):
        if not np.isfinite(reference_density) or reference_density<=0 or not isinstance(quadrature,int) or quadrature<2:
            raise ValueError('Positive fixed reference density and at least two Gauss points required')
        self.source=source;self.density=float(reference_density);self.quadrature=quadrature
        self.number=np.r_[np.arange(source.network.count),source.network.electrons]
        self.charge=np.r_[source.network.charge,-1.]
        self.frame=np.asarray(getattr(source,'reference_velocity',np.zeros(3)))

    def from_material(self,model,state,species):
        """Gather physical moments without enforcing neutrality or common support."""
        species=np.asarray(species)
        if (species.shape!=(model.species,) or not np.issubdtype(species.dtype,np.integer) or
                np.any(species<0) or np.any(species>=len(self.number))):
            raise ValueError('Atomic species index for every material support block required')
        if not np.array_equal(model.mass,self.source.mass[species]) or not np.array_equal(model.charge,self.charge[species]):
            raise ValueError('Material and reaction species must have identical physical masses and charges')
        model.validate(state)
        edges=np.unique(np.r_[model.bounds,state.faces.ravel()]);nodes,weights=np.polynomial.legendre.leggauss(self.quadrature)
        points=(edges[:-1,None]+edges[1:,None])/2+np.diff(edges)[:,None]*nodes/2
        measure=np.diff(edges)[:,None]*weights/2;sample=model.sample(state,points.ravel())
        rows=[];network=self.source.network
        for k in range(points.size):
            atomic=np.zeros(self.source.base_size);temperatures=np.zeros(network.count+1);velocity=np.zeros((network.count+1,3))
            active=sample['number'][:,k]>0;ids=species[active]
            if len(set(ids))!=len(ids):raise ValueError('Overlapping blocks of the same fluid require an explicit mixing closure')
            atomic[self.number[ids]]=sample['number'][active,k]/self.density
            temperatures[ids]=sample['temperature_ev'][active,k];velocity[ids]=sample['velocity'][active,k]
            atomic[network.thermal]=1.5*atomic[network.electrons]*temperatures[-1]
            rows.append(self.source.initial(atomic,temperatures[:-1],velocity))
        return ReactionOverlapState(edges,points,measure,np.asarray(rows).reshape(points.shape+(self.source.size,)),model.left_electric)

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite reaction step required')
        evolved=np.array([self.source.step(row,dt,self.density) if np.any(row) else row.copy()
                          for row in state.moments.reshape(-1,self.source.size)])
        return ReactionOverlapState(state.edges,state.points,state.weights,evolved.reshape(state.moments.shape),state.left_electric)

    def fields(self,state):
        rows=state.moments.reshape(-1,self.source.size);values=[self.source.fields(row) for row in rows]
        numbers=np.stack([v['numbers'] for v in values])*self.density
        velocity=np.stack([v['velocity'] for v in values])+self.frame
        occupied=numbers>0;velocity=np.where(occupied[...,None],velocity,0.)
        temperature=np.stack([v['temperature_ev'] for v in values])
        current=EV*np.sum(numbers[...,None]*self.charge[None,:,None]*velocity,axis=1)
        return dict(number_m3=numbers,velocity_m_s=velocity,temperature_ev=temperature,current_a_m2=current,
                    charge_c_m3=EV*(numbers@self.charge))

    def totals(self,state):
        # These are linear extensive moments; no kinetic/thermal redistribution.
        f=self.fields(state);source=self.source;network=source.network
        sums=self.density*np.sum(state.weights[...,None]*state.moments,axis=(0,1));numbers=sums[self.number]
        relative_p=source.momentum_scale*sums[source.momentum].sum(axis=0);mass=source.mass@numbers
        energy=EV*(sums[source.energy].sum()+network.energy_ev@sums[:network.count]+sums[network.photons]-sums[network.heating])
        energy+=self.frame@relative_p+.5*mass*(self.frame@self.frame)
        # Reactions conserve charge locally. The initial piecewise-constant
        # charge remains so up to source roundoff; no projection enforces it.
        charge=f['charge_c_m3'].reshape(state.points.shape)
        average=np.sum(state.weights*charge,axis=1)/np.diff(state.edges)
        electric=np.r_[state.left_electric,state.left_electric+np.cumsum(average*np.diff(state.edges))/EPS0]
        field_energy=EPS0/6*np.sum(np.diff(state.edges)*(electric[:-1]**2+electric[:-1]*electric[1:]+electric[1:]**2))
        return dict(numbers_m2=numbers,nuclei_m2=float(numbers[:-1].sum()),charge_c_m2=float(EV*(numbers@self.charge)),
                    mass_kg_m2=float(mass),momentum_kg_per_m_s=relative_p+mass*self.frame,
                    matter_binding_photon_energy_j_m2=float(energy),electric_energy_j_m2=float(field_energy),
                    energy_j_m2=float(energy+field_energy),photon_energy_j_m2=float(EV*sums[network.photons]),
                    electric_edges_v_m=electric)
