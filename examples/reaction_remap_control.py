"""Full20ns zero-source control isolating repeated moving-mesh reaction return."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import atomic_mass
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_vector_moments import VectorizedAtomicMoments
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination
from plasma_bench.precise_geometry import quadrature
from plasma_bench.precise_material import PreciseMaterialSubcells, PreciseMaterialSet
from plasma_bench.material_species_set import MaterialSetState
from plasma_bench.reactive_material import MaterialReservoirs
from plasma_bench.subcycled_material import SubcycledReactiveMaterial
from .precise_material import save_states
from .implicit_plasma_records import provenance, freeze


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh output required')
    ion = ImpactIonization.from_file('data/atomic/argon/argon.json')
    rr = RadiativeRecombination.from_file('data/atomic/recombination/argon_rr.json')
    network = AtomicNetwork(ion, rr, ionization_enabled=False, recombination_enabled=False)
    source = VectorizedAtomicMoments(network, 39.948*atomic_mass, [0., 0., 0.])
    identity = provenance(['examples.reaction_remap_control'])
    out.mkdir(parents=True)
    freeze(out, identity)
    report = dict(terminal=False, source_sha256=identity['source_sha256'], cases=[],
        scope='Original20ns argon initial states with DI/RR disabled in every case. '
              'Compare direct transport with exact-zero source gather/return. '
              'Isolates numerical transfer; not a reacting experiment or source-time convergence.',
        inputs=dict(ionization_sha256=ion.input_sha256, recombination_sha256=rr.input_sha256))
    began = time.perf_counter()
    reference = None
    for interval_ns in (0., 10., 5.):
        faces = -.01*np.cos(np.linspace(0, np.pi, 33))[None, :]
        ids = [network.ground[0], network.ground[1], network.count]
        models, blocks = [], []
        for index, temperature in zip(ids, [.1, .2, 2.]):
            m = PreciseMaterialSubcells(1e14*np.diff(faces)[0], source.mass[index],
                np.r_[network.charge, -1.][index], (-.2, .2), faces[0], 32,
                maximum_plasma_phase=.04, numerical_viscosity=(.5, 2.), viscosity_limiter='minmod')
            models.append(m)
            blocks.append(m.state(faces, temperature))
        model, state = PreciseMaterialSet(models, ids), MaterialSetState(tuple(blocks))
        ledger = MaterialReservoirs()
        solver = SubcycledReactiveMaterial(source, 1e14, velocity_cells=32,
            maximum_plasma_phase=.04, coalesce_ulps=0, transport_fraction=.5)
        initial_energy = solver.totals(model, state, ledger)['full_energy_j_m2']
        elapsed, steps, history, failure = 0., 0, [], None
        name = 'direct' if not interval_ns else f'remap_{interval_ns:g}ns'
        print(name, 'begin', flush=True)
        save_states(out/(name+'_initial.npz'), [model], [state])
        while elapsed < 20e-9:
            dt = min(interval_ns*1e-9 if interval_ns else .5*model.timestep(state), 20e-9-elapsed)
            try:
                if interval_ns:
                    model, state, ledger = solver.step(model, state, ledger, dt)
                    history.extend(dict(coupled_interval_start_s=elapsed, **r) for r in solver.last_transport)
                    steps += len(solver.last_transport)
                else:
                    state = model.step(state, dt)
                    steps += 1
                    history.append(dict(time_s=elapsed+dt, dt_s=dt,
                        energy_j_m2=solver.totals(model, state, ledger)['full_energy_j_m2']))
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                failure = str(error)
                break
            elapsed = min(elapsed+dt, 20e-9)
        save_states(out/(name+'_final.npz'), [model], [state])
        if history:
            np.savez_compressed(out/(name+'_transport.npz'), **{
                key: np.asarray([r[key] for r in history]) for key in history[0]})
        row = dict(name=name, source_interval_ns=interval_ns, completed_time_s=elapsed,
            requested_time_s=20e-9, success=failure is None and elapsed == 20e-9,
            failure=failure, transport_steps=steps, occupied_density_cells=sum(m.cells for m in model.models),
            relative_energy_error=(solver.totals(model, state, ledger)['full_energy_j_m2']-initial_energy)/initial_energy)
        if row['success'] and not interval_ns:
            reference = model, state
        elif row['success'] and reference is not None:
            original, original_state = reference
            edges = np.unique(np.r_[model.bounds, state.faces, original_state.faces])
            points, weights = quadrature(edges, 2)
            a, b = model.sample(state, points.ravel()), original.sample(original_state, points.ravel())
            for key in ('electric', 'current'):
                weight = weights.ravel() if key == 'electric' else weights.ravel()[:, None]
                row[key+'_relative_l2_from_direct'] = float(np.sqrt(np.sum(weight*(a[key]-b[key])**2)/np.sum(weight*b[key]**2)))
        report['cases'].append(row)
        report['elapsed_seconds'] = time.perf_counter()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
        print(row, flush=True)
    if identity != provenance(['examples.reaction_remap_control']):
        raise RuntimeError('Source changed during zero-source control')
    report['terminal'] = True
    report['success'] = all(r['success'] for r in report['cases'])
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    args = parser.parse_args()
    result = run(args.destination)
    raise SystemExit(0 if result['success'] else 1)
