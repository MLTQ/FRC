"""Conservative numerical thermal mixing, separate from physical heat transport."""
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import expm_multiply
from .material_carrier_current import MaterialCurrentState
from .material_coil_current import MaterialCoilCurrent
from .material_split_current import CentralEntropy


def diffuse_energy(mass, tail, head, conductance, internal, duration):
    """Evolve M dh/dt=-I.T diag(g) I h without changing total thermal energy."""
    mass = np.asarray(mass, float); internal = np.asarray(internal, float)
    conductance = np.asarray(conductance, float)
    tail = np.asarray(tail, int); head = np.asarray(head, int)
    if (mass.ndim != 1 or internal.shape != mass.shape or
            conductance.shape != tail.shape or tail.shape != head.shape or
            not np.isfinite(mass).all() or not np.isfinite(internal).all() or
            not np.isfinite(conductance).all() or not np.isfinite(duration) or
            np.any(mass <= 0) or np.any(internal <= 0) or
            np.any(conductance < 0) or duration < 0):
        raise ValueError('Positive finite particles/energy and nonnegative mixing/time required')
    if duration == 0 or not np.any(conductance):
        return internal.copy()
    rows = np.concatenate((tail, head, tail, head))
    cols = np.concatenate((tail, head, head, tail))
    values = np.concatenate((conductance, conductance, -conductance, -conductance))
    laplacian = csr_matrix((values, (rows, cols)), shape=(len(mass), len(mass)))
    generator = -duration*diags(1/mass)@laplacian
    evolved = expm_multiply(generator, internal, traceA=float(generator.diagonal().sum()))
    if not np.isfinite(evolved).all() or np.any(evolved <= 0):
        raise RuntimeError('Thermal mixing lost positive finite energy; no clipping or repair')
    return evolved


class EnergyDiffusionMaterialCoilCurrent(MaterialCoilCurrent):
    """Central physical entropy/pressure RK4 with conservative mixing halves."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.species = CentralEntropy(self.fluid, self.space)
        self.split_diagnostics = {}
        self.mixing_entropy = np.zeros(2)
        self.mixing_redistribution = np.zeros(2)
        self.mixing_residual = np.zeros(2)

    def diffusion_step(self, state, duration):
        species = self.species; geometry = self.geometry(state.fluid)
        edges = np.flatnonzero(species.edge_scale > 0)
        tail, head = species.tail[edges], species.head[edges]
        flux = -species.edge_scale[edges]*np.asarray(geometry['d'].T@state.magnetic)[edges]/E
        nodes = species.state_nodes(state.fluid.entropy)
        internal = geometry['fields']['specific_internal'].reshape(-1, 2)
        h0 = np.column_stack([species.project(h) for h in internal.T])
        h1 = np.column_stack([diffuse_energy(species.node_particles, tail, head,
            .5*abs(fraction*flux), h0[:, k], duration)
            for k, fraction in enumerate(species.fractions)])
        increment = np.log(h1/h0)
        new_nodes = nodes+increment
        entropy = (species.sample@new_nodes).reshape(state.fluid.entropy.shape)
        result = MaterialCurrentState(replace(state.fluid, entropy=entropy), state.magnetic.copy())
        # Observe the energy represented by the actual logarithmic state, not
        # a repaired target. Geometry and every other state are fixed here.
        represented = h0*np.exp(new_nodes-nodes)
        change = species.node_particles[:, None]*(represented-h0)
        entropy_change = np.sum(species.node_particles[:, None]*(new_nodes-nodes)
                                /(self.fluid.gammas-1), axis=0)
        return result, dict(thermal_energy_residual_j=change.sum(axis=0),
            entropy_production_over_kB=entropy_change,
            node_energy_redistributed_j=.5*np.sum(abs(change), axis=0))

    def step(self, state, dt, time=0.):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite step required')
        first, a = self.diffusion_step(state, .5*dt)
        middle = super().step(first, dt, time)
        result, b = self.diffusion_step(middle, .5*dt)
        self.response(result, rate=False, time=time+dt)
        combined = {key: a[key]+b[key] for key in a}
        self.mixing_entropy += combined['entropy_production_over_kB']
        self.mixing_redistribution += combined['node_energy_redistributed_j']
        self.mixing_residual += combined['thermal_energy_residual_j']
        self.split_diagnostics = {key: value.tolist() for key, value in combined.items()}
        self.split_diagnostics.update(
            cumulative_entropy_production_over_kB=self.mixing_entropy.tolist(),
            cumulative_node_energy_redistributed_j=self.mixing_redistribution.tolist(),
            cumulative_thermal_energy_residual_j=self.mixing_residual.tolist())
        return result
