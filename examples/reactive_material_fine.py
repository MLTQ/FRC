"""Moving argon with retained density resolution and measured computational cost."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import sys
import time
import numpy as np
from scipy.constants import atomic_mass
from plasma_bench.atomic_vector_moments import VectorizedAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination
from plasma_bench.material_subcells import MaterialSubcells
from plasma_bench.material_species_set import MaterialSpeciesSet, MaterialSetState
from plasma_bench.reactive_material import ReactiveMaterial, MaterialReservoirs
from .reactive_material import write_states
from .implicit_plasma_records import provenance, freeze


def run(destination, cells=32):
    out = Path(destination)
    if out.exists() or cells < 2:
        raise ValueError('Fresh output and at least two velocity cells required')
    ion = ImpactIonization.from_file('data/atomic/argon/argon.json')
    rr = RadiativeRecombination.from_file('data/atomic/recombination/argon_rr.json')
    network = AtomicNetwork(ion, rr)
    source = VectorizedAtomicMoments(network, 39.948*atomic_mass, [0., 0., 0.])
    identity = provenance(['examples.reactive_material_fine'])
    out.mkdir(parents=True)
    freeze(out, identity)
    began = time.perf_counter()
    ids = [network.ground[0], network.ground[1], network.count]
    faces = -.01*np.cos(np.linspace(0, np.pi, cells+1))[None, :]
    models, blocks = [], []
    for index, temperature in zip(ids, [.1, .2, 2.]):
        m = MaterialSubcells(1e14*np.diff(faces)[0], source.mass[index],
            np.r_[network.charge, -1.][index], (-.2, .2), faces[0], cells,
            maximum_plasma_phase=.04*32/cells, numerical_viscosity=(.5, 2.), viscosity_limiter='minmod')
        models.append(m)
        blocks.append(m.state(faces, temperature))
    model = MaterialSpeciesSet(models, ids)
    state = MaterialSetState(tuple(blocks))
    ledger = MaterialReservoirs()
    solver = ReactiveMaterial(source, 1e14, velocity_cells=cells,
        maximum_plasma_phase=.04*32/cells, coalesce_ulps=8, minimum_density_width=0.)
    duration = 20e-9
    report = dict(terminal=False, source_sha256=identity['source_sha256'], velocity_cells=cells,
        requested_time_s=duration, completed_time_s=0., minimum_density_width_m=0.,
        roundoff_coalesce_ulps=8, maximum_observation_seconds=600., maximum_rss_bytes=2*2**30,
        minimum_free_disk_bytes=256*2**20, steps=[], failure=None,
        inputs=dict(ionization_sha256=ion.input_sha256, recombination_sha256=rr.input_sha256),
        scope='Full original20ns argon target with unchanged physical masses, ground DI/RR, '
              'fine density and separate P1 mechanics. No physical-width coarsening or fixed cell cap. '
              'Resource-limited prefixes are failures, not evidence that coarsening is necessary.')
    initial_energy = solver.totals(model, state, ledger)['full_energy_j_m2']
    def save():
        report['elapsed_seconds'] = time.perf_counter()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    def record(index, elapsed, phases):
        totals = solver.totals(model, state, ledger)
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
        arrays = [a for m, s in zip(model.models, state.blocks) for a in
                  (s.faces, s.velocity, s.thermal, m.relative_numbers, m.reference_faces, m.velocity_indices)]
        estimate = sum(a.nbytes for a in arrays)+128*1024
        free = shutil.disk_usage(out).free
        row = dict(step=index, time_s=elapsed, occupied_density_cells=sum(m.cells for m in model.models),
            peak_rss_bytes=rss, free_disk_bytes=free, estimated_uncompressed_snapshot_bytes=estimate,
            atomic_numbers_m2=totals['atomic_numbers_m2'].tolist(), full_energy_j_m2=totals['full_energy_j_m2'],
            relative_energy_change=(totals['full_energy_j_m2']-initial_energy)/initial_energy,
            photon_energy_j_m2=ledger.photons_j_m2, observed_projection_field_change_j_m2=ledger.projection_field_change_j_m2,
            smallest_cell_m=min(float(np.diff(s.faces).min()) for s in state.blocks), **phases)
        if free-estimate < report['minimum_free_disk_bytes']:
            row['snapshot_saved'] = False
            report['failure'] = 'Insufficient disk headroom for the accepted state snapshot; no shorter-duration success'
        else:
            write_states(out/f'step_{index:05}.npz', [model], [state])
            row['snapshot_saved'] = True
        report['steps'].append(row)
        report['completed_time_s'] = elapsed
        save()
        print('fine argon', index, 'time_s', elapsed, 'cells', row['occupied_density_cells'],
              'rss_MiB', rss/2**20, 'phase_seconds', phases, flush=True)
    record(0, 0., {})
    elapsed = 0.
    while elapsed < duration and report['failure'] is None:
        previous = report['steps'][-1]
        if time.perf_counter()-began > report['maximum_observation_seconds'] or previous['peak_rss_bytes'] > report['maximum_rss_bytes']:
            report['failure'] = 'Declared wall-time or memory observation bound; full20ns not completed'
            break
        dt = min(.5*model.timestep(state), duration-elapsed)
        if not np.isfinite(dt) or dt <= 0 or elapsed+dt == elapsed:
            report['failure'] = 'Material timestep cannot advance represented physical time'
            break
        try:
            started = time.perf_counter()
            new_model, new_state, new_ledger = solver.source_step(model, state, ledger, dt/2)
            first = time.perf_counter()
            if dt > new_model.timestep(new_state)*(1+1e-12):
                raise ValueError('Reaction changed the admissible material timestep; no hidden retry')
            new_state = new_model.step(new_state, dt)
            moved = time.perf_counter()
            new_model, new_state, new_ledger = solver.source_step(new_model, new_state, new_ledger, dt/2)
            final = time.perf_counter()
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            report['failure'] = str(error)
            report['rejected_step_s'] = dt
            break
        model, state, ledger = new_model, new_state, new_ledger
        elapsed = min(elapsed+dt, duration)
        record(len(report['steps']), elapsed, dict(source_first_s=first-started,
               material_s=moved-first, source_second_s=final-moved))
    if provenance(['examples.reactive_material_fine']) != identity:
        raise RuntimeError('Source changed during retained-resolution calculation')
    report['terminal'] = True
    report['success'] = report['failure'] is None and elapsed == duration
    save()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    parser.add_argument('--cells', type=int, default=32)
    args = parser.parse_args()
    result = run(args.destination, args.cells)
    raise SystemExit(0 if result['success'] else 1)
