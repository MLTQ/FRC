"""Argon ground-neutral elastic exchange within the charged/reaction source."""
import numpy as np
from .atomic_vector_moments import VectorizedCollisionalAtomicMoments
from .elastic_maxwellian import TabulatedElasticPair


class ElasticAtomicMoments(VectorizedCollisionalAtomicMoments):
    """One supplied ground-neutral/electron pair plus the existing charged pairs."""
    def __init__(self, network, neutral_mass, reference_velocity, *, coulomb_log, elastic_section):
        super().__init__(network, neutral_mass, reference_velocity, coulomb_log=coulomb_log)
        if elastic_section.target != 'Ar' or network.ionization.nuclear_charge != 18:
            raise ValueError('This explicit wrapper requires an argon ground-neutral table and network')
        self.neutral_index = network.ground[0]
        self.elastic_pair = TabulatedElasticPair(elastic_section, self.mass[-1], self.mass[self.neutral_index])
        self.maximum_elastic_outside_probability = 0.

    def collision_generator(self, state, density):
        matrix, diagnostic = super().collision_generator(state, density)
        f = self.fields(state); electron, neutral = self.network.count, self.neutral_index
        diagnostic['electron_ground_neutral_collision_pairs'] = 0
        if f['numbers'][electron] > 0 and f['numbers'][neutral] > 0:
            c = self.elastic_pair.coefficients(f['temperature_ev'][electron], f['temperature_ev'][neutral],
                                              np.linalg.norm(f['velocity'][electron]-f['velocity'][neutral]))
            self.add_pair_generator(matrix, f, density, np.array([electron]), np.array([neutral]),
                                    np.array([c['drag']]), np.array([c['exchange']]))
            diagnostic['electron_ground_neutral_collision_pairs'] = 1
            diagnostic['elastic_probability_outside_table'] = c['probability_outside_table']
            self.maximum_elastic_outside_probability = max(self.maximum_elastic_outside_probability,
                                                           c['probability_outside_table'])
        return matrix, diagnostic
