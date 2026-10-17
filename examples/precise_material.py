"""Full analytic expansion/decay with actual precise material geometry."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, electron_mass as ME, elementary_charge as EV, epsilon_0 as EPS0
from plasma_bench.precise_geometry import ARITHMETIC, BITS, DIGITS, coordinates, widths
from plasma_bench.precise_material import PreciseMaterialSubcells, PreciseMaterialSet
from plasma_bench.precise_reactions import PreciseReactiveMaterial
from plasma_bench.subcycled_material import SubcycledReactiveMaterial
from plasma_bench.material_species_set import MaterialSetState
from plasma_bench.reactive_material import MaterialReservoirs
from plasma_bench.atomic_vector_moments import VectorizedAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from .reactive_material_affine import WarmDecay
from .reactive_material import recombination
from .atomic_moments import ConstantIonization
from .material_subcells import subdivided_faces
from .atomic_moving import charged_expansion_time
from .implicit_plasma_records import provenance, freeze


def save_states(path, models, states):
    arrays = {}
    for k, (model, state) in enumerate(zip(models, states)):
        for b, (m, s) in enumerate(zip(model.models, state.blocks)):
            p = f'step{k}_block{b}_'
            arrays.update({p+'faces_decimal_m': s.faces.astype(str),
                p+'reference_faces_decimal_m': m.reference_faces.astype(str),
                p+'relative_velocity_m_s': s.velocity, p+'thermal_divided_by_number_scale': s.thermal,
                p+'numbers_divided_by_number_scale': m.relative_numbers,
                p+'number_scale_m2': np.asarray(m.number_scale), p+'velocity_indices': m.velocity_indices,
                p+'extend_velocity': np.asarray(m.extend_velocity)})
        arrays[f'step{k}_species_indices'] = model.species_indices
        arrays[f'step{k}_frame_velocity_m_s'] = model.reference_velocity
    np.savez_compressed(path, **arrays)
    # Reload the physical endpoint coordinates as text, never pickled objects.
    with np.load(path, allow_pickle=False) as saved:
        for b, s in enumerate(states[-1].blocks):
            if not np.array_equal(coordinates(saved[f'step{len(states)-1}_block{b}_faces_decimal_m']), s.faces):
                raise ValueError('Saved precise endpoint coordinates did not round-trip')


def run(destination, resolutions, projection_method='galerkin', extend_velocity=False,
        newborn_spacing_fraction=0., velocity_refinement_rounds=0,
        maximum_refined_velocity_cells=256, subcycle_transport=False):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh precise-material output required')
    identity = provenance(['examples.precise_material'])
    out.mkdir(parents=True)
    freeze(out, identity)
    began = time.perf_counter()
    network = AtomicNetwork(ConstantIonization(), recombination(), {0: WarmDecay()},
                            ionization_enabled=False, recombination_enabled=False)
    source = VectorizedAtomicMoments(network, 1.00794*atomic_mass, [0., 0., 0.])
    report = dict(terminal=False, source_sha256=identity['source_sha256'], geometry_decimal_digits=DIGITS,
        geometry_arithmetic=ARITHMETIC, geometry_precision_bits=BITS,
        projection_method=projection_method,
        extend_velocity=extend_velocity,
        newborn_spacing_fraction=newborn_spacing_fraction,
        velocity_refinement_rounds=velocity_refinement_rounds,
        maximum_refined_velocity_cells=maximum_refined_velocity_cells,
        subcycle_transport=bool(subcycle_transport),
        minimum_density_width_m=0., coalesce_ulps=8, cases=[],
        scope='Full warm neutral affine expansion with actual radiative transfer plus independent charged '
              'affine expansion. Mathematical component limits, not physical argon/FRC formation.')
    def save():
        report['elapsed_seconds'] = time.perf_counter()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for kind in ('decay', 'charged'):
        for cells in resolutions:
            mass, density, temperature = (source.mass[1], 1e19, .1) if kind == 'decay' else (39.948*atomic_mass-ME, 1e12, 1.)
            a = .025
            faces = subdivided_faces(cells, a)
            thermal_profile = temperature*(1-(faces[:-1]**2+faces[:-1]*faces[1:]+faces[1:]**2)/(3*a*a))
            left_field = 0. if kind == 'decay' else -EV*density*a/EPS0
            m = PreciseMaterialSubcells(density*np.diff(faces), mass, 0. if kind == 'decay' else 1.,
                (-.2, .2), faces, cells, left_electric=left_field,
                numerical_viscosity=(.5, 2.), viscosity_limiter='minmod')
            model = PreciseMaterialSet([m], [1])
            state = MaterialSetState((m.state(faces[None, :], thermal_profile),))
            ledger = MaterialReservoirs()
            solver_type = SubcycledReactiveMaterial if subcycle_transport else PreciseReactiveMaterial
            solver = solver_type(source, density, velocity_cells=cells, coalesce_ulps=8,
                projection_method=projection_method, extend_velocity=extend_velocity,
                newborn_spacing_fraction=newborn_spacing_fraction,
                velocity_refinement_rounds=velocity_refinement_rounds,
                maximum_refined_velocity_cells=maximum_refined_velocity_cells)
            duration, rate = charged_expansion_time(1.5, 2*EV*temperature/(mass*a*a),
                0. if kind == 'decay' else EV**2*density/(EPS0*mass))
            totals = (lambda: solver.totals(model, state, ledger)) if kind == 'decay' else (lambda: model.totals(state))
            times, models, states, budgets = [0.], [model], [state], [totals()]
            transport_times, transport_dts, failed_transport = [], [], []
            failure = None
            while times[-1] < duration:
                dt = min(.5*model.timestep(state), duration-times[-1])
                if not np.isfinite(dt) or dt <= 0 or times[-1]+dt == times[-1]:
                    failure = 'Actual material step cannot advance represented physical time'
                    break
                if len(times) > 10000 or sum(q.cells for q in model.models) > 100000:
                    failure = 'Declared observation budget reached; original duration not complete'
                    break
                try:
                    if kind == 'decay':
                        new_model, new_state, new_ledger = solver.step(model, state, ledger, dt)
                        if subcycle_transport:
                            transport_times.extend(times[-1]+r['time_within_transport_s'] for r in solver.last_transport)
                            transport_dts.extend(r['dt_s'] for r in solver.last_transport)
                        else:
                            transport_times.append(times[-1]+dt); transport_dts.append(dt)
                        model, state, ledger = new_model, new_state, new_ledger
                    else:
                        state = model.step(state, dt)
                        transport_times.append(times[-1]+dt); transport_dts.append(dt)
                except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError) as error:
                    failure = str(error)
                    if kind == 'decay' and subcycle_transport:
                        failed_transport = solver.last_transport.copy()
                    break
                times.append(min(times[-1]+dt, duration))
                models.append(model)
                states.append(state)
                budgets.append(totals())
                if len(times) == 2 or (len(times)-1) % 50 == 0:
                    print(kind, cells, 'step', len(times)-1, 'time', times[-1],
                          'cells', sum(q.cells for q in model.models), flush=True)
            energy = np.array([r['full_energy_j_m2' if kind == 'decay' else 'energy_j_m2'] for r in budgets])
            row = dict(name=f'{kind}_n{cells}', success=failure is None, failure=failure,
                steps=len(times)-1, requested_time_s=duration, completed_time_s=times[-1],
                actual_transport_steps=len(transport_dts), failed_transport_records=failed_transport,
                maximum_occupied_cells=max(sum(m.cells for m in q.models) for q in models),
                minimum_occupied_width_m=min(float(widths(s.faces).min()) for q in states for s in q.blocks),
                numerical_projection_thermal_change_j_m2=ledger.projection_thermal_gain_j_m2,
                projection_energy_redistribution_l1_j_m2=ledger.projection_energy_redistribution_l1_j_m2,
                maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]))
            extra = {}
            if kind == 'decay':
                number = np.array([r['atomic_numbers_m2'] for r in budgets])
                photons = np.array([r['photon_energy_j_m2'] for r in budgets])
                exact = 2*a*density*np.exp(-2e5*np.asarray(times))
                row.update(maximum_relative_population_error=float(abs(number[:, 1]-exact).max()/(2*a*density)),
                    maximum_photon_error_j_m2=float(abs(photons-3*EV*(2*a*density-exact)).max()))
                extra.update(atomic_numbers_m2=number, photon_energy_j_m2=photons)
                axes[0].plot(np.asarray(times)*1e6, number[:, 1]/(2*a*density), label=f'{cells} cells')
            if failure is None:
                row.update(maximum_front_error_m=max(float(abs(np.asarray(s.faces[0, [0, -1]], float)-1.5*np.array([-a, a])).max()) for s in state.blocks),
                    maximum_velocity_error_m_s=max(float(abs(s.velocity[0, :, 2]-rate*np.asarray(s.faces[0, q.velocity_indices], float)/1.5).max()) for q, s in zip(model.models, state.blocks)))
            save_states(out/(row['name']+'.npz'), models, states)
            np.savez_compressed(out/(row['name']+'_budgets.npz'), time_s=times, energy_j_m2=energy,
                transport_stage_time_s=transport_times, transport_dt_s=transport_dts, **extra)
            report['cases'].append(row)
            save()
            print(row, flush=True)
            axes[1 if kind == 'decay' else 2].plot(np.asarray(times)*1e6, (energy-energy[0])/energy[0],
                label=f'{cells} cells'+(' STOPPED' if failure else ''))
    neutral_duration, _ = charged_expansion_time(1.5, 2*EV*.1/(source.mass[1]*.025**2), 0.)
    t = np.linspace(0, neutral_duration, 200)
    axes[0].plot(t*1e6, np.exp(-2e5*t), 'k--', label='Analytic decay')
    for ax in axes:
        ax.set_xlabel('Time [µs]')
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    axes[0].set_ylabel('Excited population / initial')
    axes[1].set_ylabel('Decay: relative total-energy error')
    axes[2].set_ylabel('Charged: relative total-energy error')
    fig.savefig(out/'precise_expansion.png', dpi=160)
    fig.savefig(out/'precise_expansion.pdf')
    plt.close(fig)
    if identity != provenance(['examples.precise_material']):
        raise RuntimeError('Source changed during precise-geometry evolution')
    report['terminal'] = True
    report['success'] = all(r['success'] for r in report['cases'])
    save()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    parser.add_argument('--cells', nargs='+', type=int, default=[16, 32, 64])
    parser.add_argument('--projection-method', choices=('galerkin', 'bounded_moments', 'constrained_velocity', 'variance'), default='galerkin')
    parser.add_argument('--extend-velocity', action='store_true')
    parser.add_argument('--newborn-spacing-fraction', type=float, default=0.)
    parser.add_argument('--velocity-refinement-rounds', type=int, default=0)
    parser.add_argument('--maximum-refined-velocity-cells', type=int, default=256)
    parser.add_argument('--subcycle-transport', action='store_true')
    args = parser.parse_args()
    result = run(args.destination, args.cells, args.projection_method, args.extend_velocity,
                 args.newborn_spacing_fraction, args.velocity_refinement_rounds,
                 args.maximum_refined_velocity_cells, args.subcycle_transport)
    raise SystemExit(0 if result['success'] else 1)
