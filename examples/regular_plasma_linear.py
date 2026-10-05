"""Independent stationary linearization of bounded regular two-species equations."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.species_projection import SpeciesProjection
from plasma_bench.variable_compatible_current import weak_vector


class StationaryLinearPlasma:
    def __init__(self, system, state, bias):
        self.system = system; self.d = system.fluid; self.state = state
        f = self.d.fields(state)
        if np.any(state[..., 1:4] != 0) or np.any(state[..., 0] != state[0, 0, 0, 0]):
            raise ValueError('Uniform resting linearization required')
        self.n = state[0, 0, 0, 0]
        self.entropy = f['entropy'][0, 0, 0]
        if not np.all(f['entropy'] == self.entropy):
            raise ValueError('Uniform physical species entropy required')
        self.h = f['specific_internal'][0, 0, 0]
        self.model = RegularCurrent(system.mapping, state[..., 0], system.ion_mass)
        self.solver = KernelInduction(self.model, **system.solver_options)
        self.projection = SpeciesProjection(self.model)
        self.bias = np.array([0., 0., bias])

    def derivative(self, perturbation, magnetic):
        d = self.d; delta = np.asarray(perturbation, float)
        # Linear p_s = (gamma_s-1) h_s [delta C_s+(gamma_s-s_s) delta n].
        dp = (d.thermal.gammas-1)*self.h*(delta[..., 4:]+(d.thermal.gammas-self.entropy)*delta[..., 0, None])
        pressure = dp.sum(axis=-1)
        q = delta[..., 1].copy(); q[0, 0] = d.sbp.axis_trace@q[0, 1:]
        nr = -2*np.einsum('ij,ejz->eiz', d.sbp.q, d.x*q)
        qr = -np.einsum('ij,eijz->eiz', d.sbp.q, pressure[:, None, :, :]-pressure[:, :, None, :])/d.mass
        flux_n = np.zeros_like(q[:, [0, -1]])
        flux_p = pressure[:, [0, -1]].copy()
        flux_n[:-1, 1] = flux_n[1:, 0] = d.x[:-1, -1]*(q[:-1, -1]+q[1:, 0])
        flux_p[:-1, 1] = flux_p[1:, 0] = .5*(pressure[:-1, -1]+pressure[1:, 0])
        sign = np.array([-1., 1.])[None, :, None]
        nr[:, [0, -1]] -= sign*(flux_n-2*d.x[:, [0, -1]]*q[:, [0, -1]])
        qr[:, [0, -1]] -= sign*(flux_p-pressure[:, [0, -1]])/d.mass
        rate = np.zeros_like(delta); rate[..., 0] = nr/d.radial_weight
        rate[..., 1] = qr/d.radial_weight
        for offset, coefficient in enumerate(d.coefficients, 1):
            rate[..., 0] -= coefficient/d.dz*(np.roll(delta[..., 3], -offset, axis=2)-np.roll(delta[..., 3], offset, axis=2))
            rate[..., 3] -= coefficient/(d.mass*d.dz)*(np.roll(pressure, -offset, axis=2)-np.roll(pressure, offset, axis=2))
        rate[..., 4:] = self.entropy*rate[..., 0, None]
        rate[0, 0, :, 1:3] = 0.
        geo = self.system.geometry; mapping = self.system.mapping
        current = self.model.current(magnetic); j = geo.sample_electric(current)
        velocity_state = self.state.copy(); velocity_state[..., 1:4] = delta[..., 1:4]
        velocity = mapping.velocity(velocity_state)
        force = np.cross(j, self.bias)
        rate += mapping.force_rate(self.state, force)
        forcing = weak_vector(geo, -np.cross(velocity, self.bias)+self.model.delta*force/(E*self.n))
        solved = self.solver.solve(forcing)
        potential = sum(a/(E*self.n)*self.projection.apply(dp[..., k])
                        for k, a in enumerate((-self.model.beta, self.model.alpha)))
        return rate, solved.magnetic, solved.electric+geo.gradient(potential)
