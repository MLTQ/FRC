"""Radially streamed independent physical rate and representation audit."""
import hashlib
import json
import math
import numpy as np
from .current_transport_sampling import TransportSampler
from .continuum_wall_probe import physical
from .radial_term_sampling import radial_slice


NAMES = ('density', 'radial_velocity', 'azimuthal_velocity', 'axial_velocity',
         'electron_entropy', 'ion_entropy', 'radial_B', 'azimuthal_B', 'axial_B')


class NormAccumulator:
    """Keep weighted squared integrals by radial row, including both end bands."""
    def __init__(self):
        self.rows = {region: [[], []] for region in ('volume', 'axis_element', 'wall_element')}
        self.maximum = 0.

    def add(self, sample, value, target):
        error = value-target
        vector = error.ndim == sample.weight.ndim+1
        square = lambda x: np.sum(x*x, axis=-1) if vector else x*x
        magnitude = np.linalg.norm(error, axis=-1) if vector else abs(error)
        self.maximum = max(self.maximum, float(np.max(magnitude)))
        numerator = np.sum(sample.weight*square(error), axis=1)
        denominator = np.sum(sample.weight*square(target), axis=1)
        masks = dict(volume=np.ones(len(sample.r), dtype=bool),
            axis_element=sample.r <= sample.d.radius/sample.d.elements,
            wall_element=sample.r >= sample.d.radius*(1-1/sample.d.elements))
        for region, mask in masks.items():
            self.rows[region][0].extend(numerator[mask].tolist())
            self.rows[region][1].extend(denominator[mask].tolist())

    def result(self):
        values = {key: float(np.sqrt(math.fsum(a)/max(math.fsum(b), np.finfo(float).tiny)))
                  for key, (a, b) in self.rows.items()}
        values['relative_L2'] = values.pop('volume')
        return dict(relative_L2=values.pop('relative_L2'), max_absolute=self.maximum, **values)


class FieldFingerprints:
    """Hash ordered row chunks using the original whole-array shape header."""
    def __init__(self, total_rows):
        self.total_rows = total_rows
        self.hashes = {}
        self.counts = {}

    def add(self, **fields):
        for name, value in fields.items():
            value = np.ascontiguousarray(value)
            if name not in self.hashes:
                shape = (self.total_rows, *value.shape[1:])
                header = json.dumps(dict(shape=shape, dtype=value.dtype.str), sort_keys=True).encode()
                self.hashes[name] = hashlib.sha256(header)
                self.counts[name] = 0
            self.hashes[name].update(value.tobytes())
            self.counts[name] += len(value)

    def result(self):
        if any(count != self.total_rows for count in self.counts.values()):
            raise RuntimeError('Incomplete physical field fingerprint')
        return {name: value.hexdigest() for name, value in self.hashes.items()}


def measure(mode, d, system, u, du, db, epsilon, rates, order, chunk_cells=2):
    if type(chunk_cells) is not int or chunk_cells < 1:
        raise ValueError('Positive integer radial cell batch required')
    sample = TransportSampler(d, system.geometry.mesh, order)
    errors = {name: NormAccumulator() for name in (*NAMES, 'electric', 'current_rate')}
    initial_errors = {name: NormAccumulator() for name in NAMES}
    prints = FieldFingerprints(len(sample.r))
    minus = u-epsilon*du; plus = u+epsilon*du
    electric_coefficients = (rates[1].electric-rates[0].electric)/(2*epsilon)
    current_coefficients = (rates[1].current-rates[0].current)/(2*epsilon)
    for start in range(0, len(sample.r), chunk_cells*order):
        s = radial_slice(sample, start, start+chunk_cells*order)
        target = mode.fields(s.R, s.Z)
        differential = (physical(s, plus, rates[1])-physical(s, minus, rates[0]))/(2*epsilon)
        electric = s.electric(electric_coefficients); current = s.electric(current_coefficients)
        expected = mode.growth*target['physical']
        for k, name in enumerate(NAMES):
            errors[name].add(s, differential[..., k], expected[..., k])
        errors['electric'].add(s, electric, target['electric'])
        errors['current_rate'].add(s, current, mode.growth*target['current'])
        n, _, _, dn, dv = s.fluid(u, du)
        entropy = np.stack([s.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
        ds = np.stack([(s.scalar(du[..., k])-entropy[..., k-4]*dn)/n for k in (4, 5)], axis=-1)
        initial = np.concatenate((dn[..., None], dv, ds, s.magnetic(db)), axis=-1)
        for k, name in enumerate(NAMES):
            initial_errors[name].add(s, initial[..., k], target['physical'][..., k])
        prints.add(physical_rate=differential, electric=electric, current_rate=current,
            target_physical=expected, target_electric=target['electric'],
            target_current_rate=mode.growth*target['current'], weight=s.weight, initial_physical=initial)
    return dict(order=order, rate_errors={k: v.result() for k, v in errors.items()},
        initial_errors={k: v.result() for k, v in initial_errors.items()}, field_fingerprints=prints.result())
