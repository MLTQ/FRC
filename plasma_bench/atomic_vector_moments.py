"""Array assembly of the existing finite-mass atomic moment generator."""
import numpy as np
from scipy.constants import elementary_charge as EV, electron_mass as ME
from .atomic_frame import FixedFrameAtomicMoments
from .atomic_collisions import CollisionalAtomicMoments


class VectorizedAtomicMoments(FixedFrameAtomicMoments):
    """Same source equations and exponential midpoint; fewer Python operations."""

    def generator(self, state, density, heat=0.):
        n = self.network
        f = self.fields(state)
        base, diagnostic = n.generator(f['atomic_state'], density, heat)
        matrix = np.zeros((self.size+1, self.size+1))
        matrix[:self.base_size, :self.base_size] = base[:-1, :-1]
        matrix[:self.base_size, -1] = base[:-1, -1]
        ve = f['velocity'][-1]
        ve2 = float(ve@ve)
        captures = (n.recombination.coefficients(f['temperature_ev'][-1]*EV)[1]/EV
                    if f['numbers'][-1] else np.zeros(n.ionization.nuclear_charge))
        transitions = base[:n.count, :n.count].copy()
        np.fill_diagonal(transitions, 0.)
        if np.any(transitions < 0):
            raise ValueError('Atomic transition rates must be nonnegative')
        destination, origin = np.nonzero(transitions)
        rate = transitions[destination, origin]
        change = n.charge[destination]-n.charge[origin]
        if np.any((change < -1) | (change > 1)):
            raise ValueError('Only adjacent charge changes have specified moment closure')
        before, after = self.mass[origin], self.mass[destination]
        po, pd = self.momentum[origin], self.momentum[destination]
        eo, ed = self.heavy_energy[origin], self.heavy_energy[destination]

        # Repeated origins are intentional: all loss channels must accumulate.
        np.add.at(matrix, (po, po), -rate[:, None])
        np.add.at(matrix, (eo, eo), -rate)
        transfer = change >= 0
        r = rate[transfer]*after[transfer]/before[transfer]
        np.add.at(matrix, (pd[transfer], po[transfer]), r[:, None])
        np.add.at(matrix, (ed[transfer], eo[transfer]), r)
        ion = change == 1
        r = rate[ion]*ME/before[ion]
        np.add.at(matrix, (self.momentum[-1], po[ion]), r[:, None])
        np.add.at(matrix, (n.thermal, eo[ion]), r)

        capture = change == -1
        r, old, new = rate[capture], before[capture], after[capture]
        o, p, q = origin[capture], po[capture], pd[capture]
        e, g = eo[capture], ed[capture]
        captured_energy = captures[n.charge[o]-1]
        electron_p = r[:, None]*ME*ve/self.momentum_scale
        np.add.at(matrix, (q, p), r[:, None])
        np.add.at(matrix, (q, o[:, None]), electron_p)
        np.add.at(matrix, (self.momentum[-1], o[:, None]), -electron_p)
        np.add.at(matrix, (g, e), r*old/new)
        cross = (r*ME/new*self.momentum_scale/EV)[:, None]*ve
        np.add.at(matrix, (g[:, None], p), cross)
        np.add.at(matrix, (g, o), r*(.5*ME**2/new*ve2/EV+ME/new*captured_energy))
        np.add.at(matrix, (n.thermal, o), -r*.5*ME*ve2/EV)
        np.add.at(matrix, (n.photons, e), r*ME/new)
        np.add.at(matrix, (n.photons, p), -cross)
        np.add.at(matrix, (n.photons, o), r*(.5*old*ME/new*ve2/EV-ME/new*captured_energy))

        occupied = f['numbers'][:-1] > 0
        if f['numbers'][-1]:
            speed = np.sqrt(2*EV*f['temperature_ev'][-1]/ME)
            relative = np.linalg.norm(f['velocity'][:-1][occupied]-ve, axis=1)/speed
            thermal = np.sqrt(ME*f['temperature_ev'][:-1][occupied]/
                              (self.mass[:-1][occupied]*f['temperature_ev'][-1]))
            diagnostic.update(maximum_drift_over_electron_thermal=float(relative.max(initial=0.)),
                              maximum_heavy_thermal_over_electron_thermal=float(thermal.max(initial=0.)))
        else:
            diagnostic.update(maximum_drift_over_electron_thermal=None,
                              maximum_heavy_thermal_over_electron_thermal=None)
        return matrix, diagnostic


class VectorizedCollisionalAtomicMoments(CollisionalAtomicMoments):
    """Existing paired collision split with array-assembled reaction moments."""
    def __init__(self, network, neutral_mass, reference_velocity, *, coulomb_log):
        super().__init__(network, neutral_mass, reference_velocity, coulomb_log=coulomb_log)
        self.reactions = VectorizedAtomicMoments(network, neutral_mass, reference_velocity)
