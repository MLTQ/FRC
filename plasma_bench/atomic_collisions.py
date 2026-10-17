"""Joint atomic reactions and paired charged-species momentum/energy exchange."""
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from scipy.special import gammainc
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from .atomic_frame import FixedFrameAtomicMoments


def pair_coefficients(ta, tb, va, vb, ma, mb, za, zb, coulomb_log):
    """Density-free drag K/(na*nb) and COM thermal exchange H/(na*nb).

    Vectorized isotropic Landau moments, including the one-cold-species limit.
    Both cold with zero relative velocity is singular and rejects explicitly.
    """
    variance = EV*(ta/ma+tb/mb); speed2 = np.sum((va-vb)**2, axis=1)
    if np.any(variance < 0) or np.any((variance == 0)&(speed2 == 0)):
        raise ValueError('Coulomb moments require finite relative thermal or drift motion')
    x2 = np.divide(speed2, 2*variance, out=np.full_like(variance, np.inf), where=variance>0)
    small = x2 < .05**2
    drag = np.empty_like(variance)
    s = variance[small]; x = x2[small]
    drag[small] = np.sqrt(2/np.pi)/(3*s**1.5)*(1+x*(-3/5+x*(3/14+x*(-1/18+x/88))))
    drag[~small] = gammainc(1.5, x2[~small])/speed2[~small]**1.5
    total = ma+mb; reduced = ma*mb/total
    coupling = 4*np.pi*(za*zb*EV**2/(4*np.pi*EPS0))**2*coulomb_log/reduced
    exchange = np.zeros_like(variance); warm = variance>0
    exchange[warm] = coupling[warm]*np.sqrt(2/np.pi)*np.exp(-x2[warm])/(total[warm]*variance[warm]**1.5)
    drag *= coupling
    if not np.isfinite(drag).all() or not np.isfinite(exchange).all():
        raise ValueError('Nonfinite physical Coulomb coefficients; no cutoff or floor')
    return drag, exchange


class CollisionalAtomicMoments(FixedFrameAtomicMoments):
    """Symmetric reaction/collision source with all charged pairs retained."""
    def __init__(self, network, neutral_mass, reference_velocity, *, coulomb_log):
        super().__init__(network, neutral_mass, reference_velocity)
        if not np.isscalar(coulomb_log) or not np.isfinite(coulomb_log) or coulomb_log<=0:
            raise ValueError('An explicit positive finite Coulomb logarithm is required')
        self.coulomb_log = float(coulomb_log)
        self.species_charge = np.r_[network.charge, -1]
        self.reactions = FixedFrameAtomicMoments(network, neutral_mass, reference_velocity)

    def collision_generator(self, state, density):
        if not np.isfinite(density) or density<=0:
            raise ValueError('Positive finite physical nucleus density required')
        matrix = np.zeros((self.size+1, self.size+1)); diagnostic = {}
        f = self.fields(state); number = f['numbers']; velocity = f['velocity']; temperature = f['temperature_ev']
        charged = np.flatnonzero((self.species_charge!=0)&(number>0))
        i, j = np.triu_indices(len(charged), 1); a, b = charged[i], charged[j]
        diagnostic['charged_collision_pairs'] = len(a)
        diagnostic['cold_charged_populations'] = int(np.count_nonzero(temperature[charged]==0))
        if not len(a):return matrix, diagnostic
        ma, mb = self.mass[a], self.mass[b]; va, vb = velocity[a], velocity[b]
        k, h = pair_coefficients(temperature[a], temperature[b], va, vb, ma, mb,
                                self.species_charge[a], self.species_charge[b], self.coulomb_log)
        self.add_pair_generator(matrix, f, density, a, b, k, h)
        return matrix, diagnostic

    def add_pair_generator(self, matrix, fields, density, a, b, k, h):
        """Accumulate conservative pair rows for supplied density-free force/heat."""
        number, velocity = fields['numbers'], fields['velocity']
        ma, mb = self.mass[a], self.mass[b]; va, vb = velocity[a], velocity[b]
        # Rates depend on the partner density, never divide by a tiny own mass density.
        nua, nub = k*density*number[b]/ma, k*density*number[a]/mb
        la, lb = (2/3)*h*density*number[b], (2/3)*h*density*number[a]
        pa, pb = self.momentum[a], self.momentum[b]; ea, eb = self.energy[a], self.energy[b]
        def momentum_add(rows, cols, values):
            np.add.at(matrix, (rows.ravel(), cols.ravel()), np.repeat(values, 3))
        momentum_add(pa, pa, -nua); momentum_add(pa, pb, nub)
        momentum_add(pb, pa, nua); momentum_add(pb, pb, -nub)
        np.add.at(matrix, (ea, ea), -la); np.add.at(matrix, (ea, eb), lb)
        np.add.at(matrix, (eb, ea), la); np.add.at(matrix, (eb, eb), -lb)
        centre = (ma[:, None]*va+mb[:, None]*vb)/(ma+mb)[:, None]
        qa = self.momentum_scale/EV*(.5*la[:, None]*va-nua[:, None]*centre)
        qb = self.momentum_scale/EV*(-.5*lb[:, None]*vb+nub[:, None]*centre)
        for row, sign in ((ea, 1.), (eb, -1.)):
            np.add.at(matrix, (np.repeat(row, 3), pa.ravel()), sign*qa.ravel())
            np.add.at(matrix, (np.repeat(row, 3), pb.ravel()), sign*qb.ravel())

    def generator(self, state, density, heat=0.):
        matrix, diagnostic = self.reactions.generator(state, density, heat)
        collisions, extra = self.collision_generator(state, density)
        diagnostic.update(extra)
        return matrix+collisions, diagnostic

    def collision_step(self, state, dt, density):
        def advance(matrix, initial, duration):
            operator = csr_matrix(duration*matrix)
            return expm_multiply(operator, np.r_[initial, 1.], traceA=operator.diagonal().sum())[:-1]
        first, _ = self.collision_generator(state, density)
        midpoint = advance(first, state, dt/2)
        second, _ = self.collision_generator(midpoint, density)
        result = advance(second, state, dt)
        self.fields(result)
        return result

    def step(self, state, dt, density, *, initial_heat=0., midpoint_heat=0.):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite step required')
        first = self.collision_step(state, dt/2, density)
        reacted = self.reactions.step(first, dt, density, initial_heat=initial_heat, midpoint_heat=midpoint_heat)
        return self.collision_step(reacted, dt/2, density)

    def observe(self, state, density):
        row = super().observe(state, density); f = self.fields(state)
        row['cold_charged_populations'] = int(np.count_nonzero((self.species_charge!=0)&(f['numbers']>0)&(f['temperature_ev']==0)))
        for label, mask in [('ion', self.species_charge>0), ('neutral', self.species_charge==0)]:
            total = f['numbers'][mask].sum()
            row[f'mean_{label}_temperature_ev'] = float(f['numbers'][mask]@f['temperature_ev'][mask]/total) if total else 0.
        return row
