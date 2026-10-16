"""One full analytic moving/decaying history for explicit source subcycling."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as EV
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_vector_moments import VectorizedAtomicMoments
from plasma_bench.precise_geometry import widths
from plasma_bench.precise_material import PreciseMaterialSubcells, PreciseMaterialSet
from plasma_bench.material_species_set import MaterialSetState
from plasma_bench.reactive_material import MaterialReservoirs
from plasma_bench.subcycled_material import SubcycledReactiveMaterial
from .reactive_material_affine import WarmDecay
from .reactive_material import recombination
from .atomic_moments import ConstantIonization
from .material_subcells import subdivided_faces
from .atomic_moving import charged_expansion_time
from .precise_material import save_states
from .implicit_plasma_records import provenance, freeze


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh subcycled affine output required')
    identity = provenance(['examples.subcycled_affine'])
    out.mkdir(parents=True)
    freeze(out, identity)
    began = time.perf_counter()
    network = AtomicNetwork(ConstantIonization(), recombination(), {0: WarmDecay()},
                            ionization_enabled=False, recombination_enabled=False)
    source = VectorizedAtomicMoments(network, 1.00794*atomic_mass, [0., 0., 0.])
    density, a, temperature, cells = 1e19, .025, .1, 16
    faces = subdivided_faces(cells, a)
    T = temperature*(1-(faces[:-1]**2+faces[:-1]*faces[1:]+faces[1:]**2)/(3*a*a))
    m = PreciseMaterialSubcells(density*np.diff(faces), source.mass[1], 0., (-.2, .2),
        faces, cells, numerical_viscosity=(.5, 2.), viscosity_limiter='minmod')
    model = PreciseMaterialSet([m], [1])
    state = MaterialSetState((m.state(faces[None, :], T),))
    ledger = MaterialReservoirs()
    solver = SubcycledReactiveMaterial(source, density, velocity_cells=cells, coalesce_ulps=8)
    duration, rate = charged_expansion_time(1.5, 2*EV*temperature/(source.mass[1]*a*a), 0.)
    times, models, states, budgets, internal = [0.], [model], [state], [solver.totals(model, state, ledger)], []
    failure = None
    for k in range(20):
        try:
            model, state, ledger = solver.step(model, state, ledger, duration/20)
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError) as error:
            failure = str(error)
            break
        internal.extend(dict(source_interval=k+1, **row) for row in solver.last_transport)
        times.append((k+1)*duration/20)
        models.append(model)
        states.append(state)
        budgets.append(solver.totals(model, state, ledger))
    energy = np.array([r['full_energy_j_m2'] for r in budgets])
    number = np.array([r['atomic_numbers_m2'] for r in budgets])
    photons = np.array([r['photon_energy_j_m2'] for r in budgets])
    exact = 2*a*density*np.exp(-2e5*np.asarray(times))
    report = dict(terminal=True, success=failure is None, failure=failure,
        source_sha256=identity['source_sha256'], requested_time_s=duration, completed_time_s=times[-1],
        source_intervals=len(times)-1, source_interval_s=duration/20, transport_steps=len(internal),
        minimum_density_width_m=0., coalesce_ulps=8,
        maximum_relative_population_error=float(abs(number[:, 1]-exact).max()/(2*a*density)),
        maximum_photon_error_j_m2=float(abs(photons-3*EV*(2*a*density-exact)).max()),
        maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),
        minimum_width_m=min(float(widths(s.faces).min()) for q in states for s in q.blocks),
        scope='Full6.0837us warm neutral affine expansion/decay with explicit20 source intervals '
              'and fine actual transport. Independent analytic motion and radiation; commuting '
              'neutral limit does not certify a charged source interval.')
    if failure is None:
        report.update(maximum_front_error_m=max(float(abs(np.asarray(s.faces[0, [0, -1]], float)-1.5*np.array([-a, a])).max()) for s in state.blocks),
            maximum_velocity_error_m_s=max(float(abs(s.velocity[0, :, 2]-rate*np.asarray(s.faces[0, q.velocity_indices], float)/1.5).max()) for q, s in zip(model.models, state.blocks)))
    save_states(out/'states.npz', models, states)
    np.savez_compressed(out/'budgets.npz', time_s=times, atomic_numbers_m2=number, energy_j_m2=energy, photon_energy_j_m2=photons)
    if internal:
        np.savez_compressed(out/'transport_stages.npz', **{key: np.array([r[key] for r in internal]) for key in internal[0]})
    if identity != provenance(['examples.subcycled_affine']):
        raise RuntimeError('Source changed during full subcycled affine trajectory')
    report['elapsed_seconds'] = time.perf_counter()-began
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(report, flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    args = parser.parse_args()
    result = run(args.destination)
    raise SystemExit(0 if result['success'] else 1)
