"""Conservative upwind species transport on material dual control volumes."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME
from scipy.sparse import csr_matrix


class MaterialDualEntropy:
    """Positive dual-volume state; explicit upwind diffusion, no clipping."""
    def __init__(self,fluid,space):
        self.fluid=fluid;self.space=space;self.particles=fluid.particles.ravel()
        corner=np.argmax(fluid.shape,axis=1);self.owner=fluid.ids[:,corner].ravel();points=len(self.particles)
        self.sample=csr_matrix((np.ones(points),(np.arange(points),self.owner)),shape=(points,fluid.nodes))
        self.node_particles=np.asarray(self.sample.T@self.particles)
        if np.any(self.node_particles<=0):raise ValueError('Positive actual particles in every dual volume required')
        gradient=space.gradient
        counts=np.count_nonzero(gradient,axis=1)
        if gradient.shape[1]!=fluid.nodes or np.any((counts!=0)&(counts!=2)):
            raise ValueError('Shared material nodes and oriented compatible edge incidence required')
        self.tail=gradient.argmin(axis=1);self.head=gradient.argmax(axis=1)
        self.edge_scale=gradient[np.arange(len(gradient)),self.head]
        self.fractions=(-fluid.ion_mass/fluid.mass,ME/fluid.mass)

    def project(self,value):
        value=np.asarray(value,float)
        if value.shape!=(len(self.particles),) or not np.isfinite(value).all():raise ValueError('Finite scalar on material quadrature required')
        return np.asarray(self.sample.T@(self.particles*value))/self.node_particles

    def state_nodes(self,entropy):
        nodes=np.column_stack([self.project(s) for s in entropy.reshape(-1,2).T])
        difference=np.max(abs(entropy.reshape(-1,2)-self.sample@nodes))
        if difference>512*np.finfo(float).eps*max(1.,float(abs(nodes).max())):
            raise ValueError('Entropy must be represented on material dual volumes; no state projection or repair')
        return nodes

    def pair(self,state,geometry,drift):
        nodes=self.state_nodes(state.entropy);internal=geometry['fields']['specific_internal'].reshape(-1,2)
        hbar=np.column_stack([self.project(h) for h in internal.T])
        particle_flux=self.edge_scale*(geometry['mass']@geometry['drift_coefficients'])
        rate=np.empty_like(nodes);forces=[];potentials=[];outflow=np.zeros_like(nodes)
        for k,(fraction,gamma) in enumerate(zip(self.fractions,self.fluid.gammas)):
            flow=fraction*particle_flux;upstream=np.where(flow>=0,self.tail,self.head)
            flux=flow*nodes[upstream,k];load=np.zeros(self.fluid.nodes)
            np.add.at(load,self.tail,-flux);np.add.at(load,self.head,flux)
            rate[:,k]=load/self.node_particles
            anchor=float(nodes[0,k]);difference=hbar[self.head,k]-hbar[self.tail,k]
            reduced=fraction/E*(nodes[upstream,k]-anchor)*self.edge_scale*difference
            forces.append(geometry['mass']@reduced)
            phi=fraction/E*(gamma*hbar[:,k]-nodes[:,k]*hbar[:,k]+anchor*hbar[:,k])
            potentials.append(self.space.gradient@phi)
            np.add.at(outflow[:,k],upstream,abs(flow))
        entropy_rate=(self.sample@rate).reshape(state.entropy.shape)
        return dict(entropy_rate=entropy_rate,reduced_forcing=sum(forces),potential_electric=sum(potentials),
            thermal_power=float(np.sum(self.node_particles[:,None]*hbar*rate)),
            entropy_rates_over_kB=np.sum(self.node_particles[:,None]*rate/(self.fluid.gammas-1),axis=0),
            maximum_transport_outflow_s=float(np.max(outflow/self.node_particles[:,None])))

    def heating_rate(self,geometry,heating):
        internal=geometry['fields']['specific_internal'].reshape(-1,2)
        denominator=self.sample.T@(self.particles[:,None]*internal)
        source=self.sample.T@(geometry['volume'][:,None]*heating)
        return (self.sample@(source/denominator)).reshape(self.fluid.particles.shape+(2,))

    def potential_gradient(self,value):return self.space.gradient@self.project(value)
