"""Time integration of the upwind pressure/current switch in a coil pulse."""
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.linalg import cho_factor, cho_solve, solve
from .material_dual_entropy import MaterialDualEntropy
from .material_coil_current import MaterialCoilCurrent
from .material_carrier_current import MaterialCurrentState
from .material_current_box import box_minimum as _box_minimum



class CentralEntropy(MaterialDualEntropy):
    """The central part of the same conservative entropy/pressure pair."""
    def pair(self, state, geometry, drift):
        nodes = self.state_nodes(state.entropy)
        internal = geometry['fields']['specific_internal'].reshape(-1, 2)
        hbar = np.column_stack([self.project(h) for h in internal.T])
        q = self.edge_scale*(geometry['mass']@geometry['drift_coefficients'])
        rate = np.zeros_like(nodes); forces = []; potentials = []
        outflow = np.zeros_like(nodes)
        for k, (fraction, gamma) in enumerate(zip(self.fractions, self.fluid.gammas)):
            flow = fraction*q; average = .5*(nodes[self.tail, k]+nodes[self.head, k])
            flux = flow*average; load = np.zeros(self.fluid.nodes)
            np.add.at(load, self.tail, -flux); np.add.at(load, self.head, flux)
            rate[:, k] = load/self.node_particles
            anchor = float(nodes[0, k]); dh = hbar[self.head, k]-hbar[self.tail, k]
            reduced = fraction/E*(average-anchor)*self.edge_scale*dh
            forces.append(geometry['mass']@reduced)
            phi = fraction/E*(gamma*hbar[:, k]-nodes[:, k]*hbar[:, k]+anchor*hbar[:, k])
            potentials.append(self.space.gradient@phi)
            np.add.at(outflow[:, k], np.where(flow >= 0, self.tail, self.head), abs(flow))
        return dict(entropy_rate=(self.sample@rate).reshape(state.entropy.shape),
            reduced_forcing=sum(forces), potential_electric=sum(potentials),
            thermal_power=float(np.sum(self.node_particles[:, None]*hbar*rate)),
            entropy_rates_over_kB=np.sum(self.node_particles[:, None]*rate/(self.fluid.gammas-1), axis=0),
            maximum_transport_outflow_s=float(np.max(outflow/self.node_particles[:, None])))


def _midpoint_current(q0, metric, k, duration, initial, previous_eta):
    """Bounded upwind trace at zero flux; retain both signs of pressure coupling."""
    coefficient = duration/(2*E)
    positive = np.flatnonzero(k > 0); negative = np.flatnonzero(k < 0)
    signs = np.sign(initial); eta = k*signs
    # The negative-k branch is nonconvex. Continue its incoming sign, then demand
    # consistency with the resulting midpoint flux; no coefficient is removed.
    for branch in range(12):
        fixed = np.zeros_like(k); fixed[negative] = k[negative]*signs[negative]
        qbar = q0-coefficient*(metric@fixed)
        eta = fixed.copy()
        if len(positive):
            bound = k[positive]; kscale = float(bound.max())
            block = metric[np.ix_(positive, positive)]
            qscale = max(float(abs(qbar).max()), coefficient*kscale*float(np.linalg.norm(block, np.inf)), np.finfo(float).tiny)
            matrix = coefficient*kscale*block/qscale; linear = qbar[positive]/qscale
            limits = bound/kscale
            guess = k[positive]*signs[positive] if previous_eta is None else previous_eta[positive]
            eta[positive] = kscale*_box_minimum(matrix, linear, limits, guess/kscale)
        q = q0-coefficient*(metric@eta)
        updated = np.sign(q)
        if np.array_equal(updated[negative], signs[negative]):
            return q, eta
        signs[negative] = updated[negative]
    raise RuntimeError('Upwind midpoint negative-pressure branch did not converge')


class SplitMaterialCoilCurrent(MaterialCoilCurrent):
    """Symmetric entropy/current diffusion halves around the physical RK4 step."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.central_species = CentralEntropy(self.fluid, self.space)
        self.split_diagnostics = {}

    def derivative(self, state, time=0.):
        upstream = self.species
        self.species = self.central_species
        try:
            return super().derivative(state, time)
        finally:
            self.species = upstream

    def diffusion_step(self, state, duration):
        species = self.species; g = self.geometry(state.fluid)
        edges = np.flatnonzero(species.edge_scale > 0)
        tail, head = species.tail[edges], species.head[edges]
        ds_operator = g['d'].toarray()[:, edges]*species.edge_scale[edges]
        operator = g['h']+self.mu/E**2*(g['d']@g['solve'](g['d'].T.toarray()))
        diagonal = np.sqrt(np.diag(operator))
        factor = cho_factor(operator/diagonal[:, None]/diagonal[None, :], lower=True)
        inverse_ds = cho_solve(factor, ds_operator/diagonal[:, None])/diagonal[:, None]
        metric = ds_operator.T@inverse_ds; metric = .5*(metric+metric.T)
        q0 = -ds_operator.T@state.magnetic/E
        old = species.state_nodes(state.fluid.entropy)
        h0 = np.column_stack([species.project(h) for h in g['fields']['specific_internal'].reshape(-1, 2).T])
        mass = species.node_particles; incidence = np.zeros((len(edges), len(mass)))
        incidence[np.arange(len(edges)), tail] = -1.
        incidence[np.arange(len(edges)), head] = 1.
        current = q0.copy(); entropy = old.copy(); scale = np.sqrt(mass); eta = None
        for iteration in range(60):
            next_entropy = np.empty_like(old)
            for j, fraction in enumerate(species.fractions):
                laplacian = incidence.T@((.5*abs(fraction*current))[:, None]*incidence)
                # Solve for the increment so a large entropy anchor cancels
                # before arithmetic, using the actual midpoint diffusion law.
                matrix = np.diag(mass)+.5*duration*laplacian
                delta = solve(matrix/scale[:, None]/scale[None, :],
                    -duration*(laplacian@(old[:, j]-old[0, j]))/scale, assume_a='pos')/scale
                next_entropy[:, j] = old[:, j]+delta
            difference = next_entropy-old
            mean_h = h0*np.divide(np.expm1(difference), difference,
                out=np.ones_like(difference), where=difference != 0)
            middle = .5*(old+next_entropy)
            k = -np.sum(abs(np.asarray(species.fractions))*(middle[head]-middle[tail])*(mean_h[head]-mean_h[tail]), axis=1)/(2*E)
            next_current, eta = _midpoint_current(q0, metric, k, duration, current, eta)
            current_error = float(abs(next_current-current).max()/max(float(abs(q0).max()), float(abs(next_current).max()), 1.))
            entropy_error = float(abs(next_entropy-entropy).max())
            current, entropy = next_current, next_entropy
            if current_error < 2e-9 and entropy_error < 2e-12:
                break
        else:
            raise RuntimeError(f'Coupled upwind diffusion did not converge: current {current_error:g}, entropy {entropy_error:g}')
        magnetic = state.magnetic+duration*(inverse_ds@eta)
        final = MaterialCurrentState(replace(state.fluid,
            entropy=(species.sample@entropy).reshape(state.fluid.entropy.shape)), magnetic)
        thermal_change = float(np.sum(mass[:, None]*h0*np.expm1(entropy-old)))
        field_change = float(.5*(magnetic-state.magnetic)@operator@(magnetic+state.magnetic))
        # This observation is never used to repair entropy, field or energy.
        self.split_diagnostics = dict(iterations=iteration+1, current_iteration_error=current_error,
            entropy_iteration_error=entropy_error, exchange_energy_residual_j=thermal_change+field_change,
            minimum_entropy_change=float((entropy.min(axis=0)-old.min(axis=0)).min()))
        if not np.isfinite(magnetic).all() or not np.isfinite(entropy).all():
            raise ValueError('Nonfinite upwind diffusion state')
        return final

    def step(self, state, dt, time=0.):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite step required')
        first = self.diffusion_step(state, .5*dt); first_diagnostics = self.split_diagnostics.copy()
        middle = super().step(first, dt, time)
        result = self.diffusion_step(middle, .5*dt)
        self.split_diagnostics = {key: (min(first_diagnostics[key], value) if key == 'minimum_entropy_change'
            else max(abs(first_diagnostics[key]), abs(value))) for key, value in self.split_diagnostics.items()}
        self.response(result, rate=False, time=time+dt)
        return result
