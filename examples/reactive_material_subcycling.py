"""Full20ns argon attempt with explicit source intervals and fine transport."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_vector_moments import VectorizedAtomicMoments
from plasma_bench.impact_ionization import ImpactIonization
from plasma_bench.radiative_recombination import RadiativeRecombination
from plasma_bench.precise_geometry import ARITHMETIC, BITS, DIGITS, widths
from plasma_bench.precise_material import PreciseMaterialSubcells, PreciseMaterialSet
from plasma_bench.material_species_set import MaterialSetState
from plasma_bench.reactive_material import MaterialReservoirs
from plasma_bench.subcycled_material import SubcycledReactiveMaterial
from .precise_material import save_states
from .implicit_plasma_records import provenance, freeze


def run(destination, source_interval_ns=10., cells=32, coalesce_ulps=0,
        projection_method='galerkin', extend_velocity=False, newborn_spacing_fraction=0.):
    out = Path(destination)
    if out.exists() or not np.isfinite(source_interval_ns) or source_interval_ns <= 0 or cells < 2:
        raise ValueError('Fresh output and positive source interval/velocity resolution required')
    if not isinstance(coalesce_ulps, int) or coalesce_ulps < 0:
        raise ValueError('Nonnegative explicit ULP cleanup required')
    ion = ImpactIonization.from_file('data/atomic/argon/argon.json')
    rr = RadiativeRecombination.from_file('data/atomic/recombination/argon_rr.json')
    network = AtomicNetwork(ion, rr)
    source = VectorizedAtomicMoments(network, 39.948*atomic_mass, [0., 0., 0.])
    identity = provenance(['examples.reactive_material_subcycling'])
    out.mkdir(parents=True)
    freeze(out, identity)
    began = time.perf_counter()
    report = dict(terminal=False, source_sha256=identity['source_sha256'], requested_time_s=20e-9,
        completed_time_s=0., source_interval_s=source_interval_ns*1e-9, transport_fraction=.5,
        velocity_cells=cells, preserve_velocity_knots=True, geometry_decimal_digits=DIGITS,
        geometry_arithmetic=ARITHMETIC, geometry_precision_bits=BITS,
        projection_method=projection_method,
        extend_velocity=bool(extend_velocity),
        newborn_spacing_fraction=newborn_spacing_fraction,
        minimum_density_width_m=0., coalesce_ulps=coalesce_ulps,
        intervals=[], failure=None, minimum_free_disk_bytes=256*2**20,
        inputs=dict(ionization_sha256=ion.input_sha256, recombination_sha256=rr.input_sha256),
        scope='Full original20ns ground-ladder argon target. Physical masses, all species, '
              'fine density cuts and precise supports retained. Source timing is explicit and '
              'unconverged; actual acoustic/plasma transport steps are recorded separately. '
              'No physical-width coarsening, shortened success, collision or formation claim.')
    def save():
        report['elapsed_seconds'] = time.perf_counter()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2,
            default=lambda value: value.item())+'\n')
    ids = [network.ground[0], network.ground[1], network.count]
    faces = -.01*np.cos(np.linspace(0, np.pi, cells+1))[None, :]
    models, blocks = [], []
    for index, temperature in zip(ids, [.1, .2, 2.]):
        m = PreciseMaterialSubcells(1e14*np.diff(faces)[0], source.mass[index],
            np.r_[network.charge, -1.][index], (-.2, .2), faces[0], cells,
            maximum_plasma_phase=.04*32/cells, numerical_viscosity=(.5, 2.), viscosity_limiter='minmod')
        models.append(m)
        blocks.append(m.state(faces, temperature))
    model = PreciseMaterialSet(models, ids)
    state = MaterialSetState(tuple(blocks))
    ledger = MaterialReservoirs()
    solver = SubcycledReactiveMaterial(source, 1e14, velocity_cells=cells,
        maximum_plasma_phase=.04*32/cells, coalesce_ulps=coalesce_ulps, transport_fraction=.5,
        projection_method=projection_method, progress_every=100, extend_velocity=extend_velocity,
        newborn_spacing_fraction=newborn_spacing_fraction)
    energy0 = solver.totals(model, state, ledger)['full_energy_j_m2']
    save_states(out/'initial.npz', [model], [state])
    save()
    elapsed, index = 0., 0
    while elapsed < report['requested_time_s']:
        if shutil.disk_usage(out).free < report['minimum_free_disk_bytes']+16*2**20:
            report['failure'] = 'Insufficient reserved disk headroom; full20ns not complete'
            break
        dt = min(report['source_interval_s'], report['requested_time_s']-elapsed)
        started = time.perf_counter()
        print('source interval', index+1, 'begin', elapsed, 'dt', dt, flush=True)
        try:
            new_model, new_state, new_ledger = solver.step(model, state, ledger, dt)
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError) as error:
            report['failure'] = str(error)
            report['failed_interval_dt_s'] = dt
            report['accepted_internal_steps_in_failed_interval'] = len(solver.last_transport)
            if solver.last_transport:
                np.savez_compressed(out/'failed_transport_stages.npz', **{
                    key: np.asarray([r[key] for r in solver.last_transport]) for key in solver.last_transport[0]})
            if solver.last_overlap is not None:
                r = solver.last_overlap
                knots = {} if solver.last_velocity_knots is None else {
                    f'velocity_knots_species{index}_decimal_m': points.astype(str)
                    for index, points in solver.last_velocity_knots.items()}
                np.savez_compressed(out/'rejected_source.npz', edges_decimal_m=r.edges.astype(str),
                    points_decimal_m=r.points.astype(str), weights_m=r.weights, moments=r.moments, **knots)
            break
        model, state, ledger = new_model, new_state, new_ledger
        elapsed = min(elapsed+dt, report['requested_time_s'])
        index += 1
        budget = solver.totals(model, state, ledger)
        save_states(out/f'interval_{index:04}.npz', [model], [state])
        np.savez_compressed(out/f'transport_{index:04}.npz', **{
            key: np.asarray([r[key] for r in solver.last_transport]) for key in solver.last_transport[0]})
        row = dict(interval=index, source_interval_s=dt, time_s=elapsed,
            elapsed_seconds=time.perf_counter()-started, transport_steps=len(solver.last_transport),
            minimum_transport_dt_s=min(r['dt_s'] for r in solver.last_transport),
            maximum_transport_dt_s=max(r['dt_s'] for r in solver.last_transport),
            occupied_density_cells=sum(m.cells for m in model.models),
            actual_velocity_cells_by_block=[m.velocity_cells for m in model.models],
            minimum_width_m=min(float(widths(s.faces).min()) for s in state.blocks),
            atomic_numbers_m2=budget['atomic_numbers_m2'].tolist(),
            full_energy_j_m2=budget['full_energy_j_m2'], relative_energy_error=(budget['full_energy_j_m2']-energy0)/energy0,
            maximum_internal_relative_energy_error=max(abs(r['energy_j_m2']-energy0)/energy0 for r in solver.last_transport),
            observed_projection_field_change_j_m2=ledger.projection_field_change_j_m2,
            numerical_projection_thermal_change_j_m2=ledger.projection_thermal_gain_j_m2,
            projection_energy_redistribution_l1_j_m2=ledger.projection_energy_redistribution_l1_j_m2,
            last_projection=solver.last_projection,
            photon_energy_j_m2=ledger.photons_j_m2)
        report['intervals'].append(row)
        report['completed_time_s'] = elapsed
        save()
        print({k: row[k] for k in ('interval', 'time_s', 'transport_steps',
            'occupied_density_cells', 'relative_energy_error')}, flush=True)
    if identity != provenance(['examples.reactive_material_subcycling']):
        raise RuntimeError('Source changed during explicit subcycled calculation')
    report['terminal'] = True
    report['success'] = report['failure'] is None and elapsed == report['requested_time_s']
    save()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    parser.add_argument('--source-step-ns', type=float, default=10.)
    parser.add_argument('--cells', type=int, default=32)
    parser.add_argument('--coalesce-ulps', type=int, default=0)
    parser.add_argument('--projection-method', choices=('galerkin', 'bounded_moments', 'constrained_velocity'), default='galerkin')
    parser.add_argument('--extend-velocity', action='store_true')
    parser.add_argument('--newborn-spacing-fraction', type=float, default=0.)
    args = parser.parse_args()
    result = run(args.destination, args.source_step_ns, args.cells, args.coalesce_ulps,
                 args.projection_method, args.extend_velocity, args.newborn_spacing_fraction)
    raise SystemExit(0 if result['success'] else 1)
