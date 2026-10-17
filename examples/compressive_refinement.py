"""Complete neutral shock or physical-mass plasma with compression refinement."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import atomic_mass, electron_mass, elementary_charge as EV
from plasma_bench.species_material import MaterialSpecies
from plasma_bench.material_species_set import MaterialSpeciesSet, MaterialSetState
from plasma_bench.compressive_refinement import refine_compressive_cells
from plasma_bench.mhd_benchmarks import sod_exact
from .implicit_plasma_records import provenance, freeze


def snapshot(path, model, state):
    arrays = dict(species_indices=model.species_indices, frame_velocity_m_s=model.reference_velocity)
    for i, (m, s) in enumerate(zip(model.models, state.blocks)):
        arrays.update({f'block{i}_faces_m': s.faces, f'block{i}_velocity_m_s': s.velocity,
                       f'block{i}_thermal_j_m2': s.thermal, f'block{i}_numbers_m2': m.numbers})
    np.savez_compressed(path, **arrays)


def run(destination, kind, cells, threshold, minimum_width, phase_cap, thermal_reconstruction='constant', refine=True):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh output required')
    identity = provenance(['examples.compressive_refinement'])
    out.mkdir(parents=True); freeze(out, identity)
    mass, gamma = 39.948*atomic_mass, 5/3
    a = .05 if kind == 'shock' else .01
    duration = .12*a/np.sqrt(EV/mass) if kind == 'shock' else 20e-9
    ids, masses, charges = ([0], [mass], [0.]) if kind == 'shock' else ([1, 19], [mass-electron_mass, electron_mass], [1., -1.])
    models, blocks = [], []
    for index, m, q in zip(ids, masses, charges):
        faces = np.linspace(-a, a, cells+1) if kind == 'shock' else -a*np.cos(np.linspace(0, np.pi, cells+1))
        centre = (faces[:-1]+faces[1:])/2
        density = 1e19*np.where(centre < 0, 1., .125) if kind == 'shock' else np.full(cells, 1e14)
        temperature = np.where(centre < 0, 1., .8) if kind == 'shock' else (.2 if index == 1 else 2.)
        component = MaterialSpecies((density*np.diff(faces))[None, :], [m], [q], (-.4, .4) if kind == 'shock' else (-.2, .2),
            maximum_plasma_phase=phase_cap, numerical_viscosity=(.5, 2.), viscosity_limiter='minmod')
        models.append(component); blocks.append(component.state(faces[None, :], temperature))
    model, state = MaterialSpeciesSet(models, ids), MaterialSetState(tuple(blocks))
    snapshot(out/'initial.npz', model, state)
    initial = model.totals(state)
    report = dict(terminal=False, source_sha256=identity['source_sha256'], kind=kind,
                  requested_time_s=duration, initial_cells=cells, threshold=threshold,
                  minimum_child_width_m=minimum_width, maximum_plasma_phase=phase_cap,
                  thermal_reconstruction=thermal_reconstruction,
                  refinement_enabled=refine,
                  scope='Complete planar component with explicit refinement; not formation or lab validation.')
    now, steps, budgets, events, failure = 0., [], [initial], [], None
    started = time.perf_counter()
    while now < duration:
        try:
            while refine:
                model, state, event = refine_compressive_cells(model, state, threshold, minimum_width,
                    species=[0] if kind == 'shock' else [19], thermal_reconstruction=thermal_reconstruction)
                if not event['changed']:
                    break
                events.append(dict(time_s=now, **event))
            dt = min(model.timestep(state), duration-now)
            if not np.isfinite(dt) or dt <= 0 or now+dt == now:
                raise ValueError('Transport clock cannot advance')
            candidate = model.step(state, dt)
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError) as error:
            failure = dict(time_s=now, message=str(error)); break
        state = candidate; now = min(now+dt, duration)
        steps.append(dt); budgets.append(model.totals(state))
        if len(steps) % 1000 == 0:
            print(kind, len(steps), now, [m.cells for m in model.models], flush=True)
    snapshot(out/'final.npz', model, state)
    energy = np.asarray([b['energy_j_m2'] for b in budgets])
    numbers = np.asarray([b['numbers_m2'] for b in budgets])
    momentum = np.asarray([b['momentum_kg_per_m_s'] for b in budgets])
    np.savez_compressed(out/'history.npz', dt_s=steps, energy_j_m2=energy, numbers_m2=numbers,
                        momentum_kg_per_m_s=momentum)
    report.update(terminal=True, success=bool(failure is None and now == duration), failure=failure,
        completed_time_s=now, steps=len(steps), elapsed_seconds=time.perf_counter()-started,
        final_cells=[m.cells for m in model.models], refinement_events=events,
        maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),
        maximum_relative_number_change=float(np.max(abs(numbers-numbers[0])/numbers[0])),
        maximum_momentum_change_kg_per_m_s=float(np.linalg.norm(momentum-momentum[0], axis=1).max()))
    if kind == 'shock' and report['success']:
        speed = np.sqrt(EV/mass); tau = now*speed/a
        _, vstar, pstar = sod_exact(np.array([0.]), tau, gamma)
        vstar, pstar = float(vstar[0]), float(pstar[0])
        shock = np.sqrt(gamma*.1/.125)*np.sqrt((gamma+1)/(2*gamma)*pstar/.1+(gamma-1)/(2*gamma))
        breaks = tau*np.array([-np.sqrt(gamma), vstar-np.sqrt(gamma)*pstar**((gamma-1)/(2*gamma)), vstar, shock])*a
        faces = state.blocks[0].faces[0]
        edges = np.unique(np.r_[-a/2, a/2, breaks, faces[(faces > -a/2) & (faces < a/2)]])
        nodes, weight = np.polynomial.legendre.leggauss(16)
        x = ((edges[:-1, None]+edges[1:, None])/2+np.diff(edges)[:, None]*nodes/2).ravel()
        w = (np.diff(edges)[:, None]*weight/2).ravel()/a
        sample = model.sample(state, x)
        actual = [sample['number'][0]/1e19, sample['velocity'][0, :, 2]/speed,
                  sample['number'][0]*sample['temperature_ev'][0]/1e19]
        errors = [float(np.sum(w*abs(v-r))) for v, r in zip(actual, sod_exact(x/a, tau, gamma))]
        report['central_normalized_l1'] = dict(zip(('density', 'velocity', 'pressure'), errors))
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    if identity != provenance(['examples.compressive_refinement']):
        raise RuntimeError('Source changed during evolution')
    print({k:v for k,v in report.items() if k != 'refinement_events'}, flush=True)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('destination'); p.add_argument('--case', choices=('shock','plasma'), required=True)
    p.add_argument('--cells', type=int, default=128); p.add_argument('--threshold', type=float, default=.01)
    p.add_argument('--minimum-child-width', type=float, required=True); p.add_argument('--phase-cap', type=float, default=.01)
    p.add_argument('--thermal-reconstruction', choices=('constant', 'limited_linear'), default='constant')
    p.add_argument('--no-refinement', action='store_true', help='Retain fixed material labels throughout the full trajectory')
    args = p.parse_args()
    result = run(args.destination, args.case, args.cells, args.threshold, args.minimum_child_width, args.phase_cap,
                 args.thermal_reconstruction, not args.no_refinement)
    raise SystemExit(0 if result['success'] else 1)
