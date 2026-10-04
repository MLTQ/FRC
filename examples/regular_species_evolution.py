"""Full-transit prescribed-current runs with independent physical error norms."""
import time
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_species_transport import FrozenSpeciesTransport
from .species_shear_reference import initialize, entropy, current_z
from .current_transport_sampling import TransportSampler
from .regular_species_study import norms


def audit(spec, arrays, order):
    ref, d, u, m, _ = initialize(spec)
    k = spec.get('species', 0)
    a = (-m.beta, m.alpha)[k]
    s = TransportSampler(d, m.mesh, order)
    n = s.scalar(u[..., 0])
    initial = s.scalar(u[..., 0]*arrays['initial_entropy'])/n
    final = s.scalar(u[..., 0]*arrays['final_entropy'])/n
    target = entropy(ref, s.R, s.Z, float(arrays['time']), a, spec.get('sign', 1))
    exact_initial = entropy(ref, s.R, s.Z, 0., a, spec.get('sign', 1))
    background = np.log(2.)+.08*(s.R/ref.radius)**2
    result = norms(final-background, target-background, s.weight*n, s.R, d)
    gamma = d.thermal.gammas[k]
    h = lambda v: E*np.exp(v)*(n/d.thermal.n0)**(gamma-1)/(gamma-1)
    h0, hf = h(initial), h(final)
    weight = s.weight*n
    j = s.electric(arrays['current'])
    expected_j = np.zeros_like(j)
    expected_j[..., 2] = current_z(ref, s.R, spec.get('sign', 1))
    return dict(signal_error=result,
                initial_signal_error=float(np.sqrt(np.sum(weight*(initial-exact_initial)**2)/np.sum(weight*(exact_initial-background)**2))),
                physical_entropy_drift_per_particle=float(np.sum(weight*(final-initial))/np.sum(weight)),
                physical_quadratic_relative_drift=float(np.sum(weight*(final**2-initial**2))/np.sum(weight*initial**2)),
                physical_internal_energy_relative_drift=float(np.sum(weight*(hf-h0))/np.sum(weight*h0)),
                prescribed_current_error=float(np.sqrt(np.sum(s.weight[..., None]*(j-expected_j)**2)/np.sum(s.weight[..., None]*expected_j**2))))


def run_case(spec):
    start = time.monotonic()
    ref, d, u, m, field = initialize(spec)
    k = spec.get('species', 0)
    current = field.current(m)
    operator = FrozenSpeciesTransport(m, current, k)
    fraction = (-m.beta, m.alpha)[k]
    speed = abs(fraction*current_z(ref, 0.)/(E*ref.number_density))
    transit = 2*ref.half_length/speed
    steps = spec.get('steps', 12*d.nz)
    dt = transit/steps
    initial = u[..., k+4]/u[..., 0]
    final = initial.copy()
    weight = d.volume*u[..., 0]
    history = []
    for i in range(steps):
        previous = final.copy()
        final = operator.step(final, dt)
        if (i+1) % max(1, steps//8) == 0 or i+1 == steps:
            history.append([dt*(i+1), np.sum(weight*(final-initial))/weight.sum(),
                            np.sum(weight*(final**2-initial**2))/np.sum(weight*initial**2)])
    arrays = dict(initial_entropy=initial, final_entropy=final, penultimate_entropy=previous,
                  current=current, time=np.asarray(transit), dt=np.asarray(dt), history=np.array(history))
    record = dict(spec=spec, status='complete', steps=steps, elapsed_seconds=time.monotonic()-start,
                  time=transit, dt=dt, grid_entropy_drift_per_particle=float(np.sum(weight*(final-initial))/weight.sum()),
                  grid_quadratic_relative_drift=float(np.sum(weight*(final**2-initial**2))/np.sum(weight*initial**2)))
    return record, arrays
