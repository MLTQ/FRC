"""Separate unchanged full spatial contributions on every saved signed probe."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.variable_compatible_current import weak_vector
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.kernel_induction import difference_curl
from .split_mode_audit import fixture
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance, freeze, sha


GROUPS = ('bulk', 'ideal_lorentz', 'hall', 'density_inertia', 'current_advection',
          'stretching_drift_stress', 'self_drift', 'electron_species', 'ion_species', 'wall_reaction')


def difference(value, target):
    norm = float(np.linalg.norm(target))
    return dict(relative_L2=float(np.linalg.norm(value-target)/norm) if norm else None,
                max_absolute=float(np.max(abs(value-target))), bitwise_equal=fingerprint(value)==fingerprint(target))


def backward(residual, magnitude):
    den = float(np.linalg.norm(magnitude)); num = float(np.linalg.norm(residual))
    return num/den if den else (0. if num == 0 else None)


def species_parts(operator, state, original):
    m = operator.model; geo = m.geometry; fields = operator.fluid.fields(state)
    result = []
    for k, (gamma, a) in enumerate(zip(operator.fluid.thermal.gammas, (-m.beta, m.alpha))):
        s = fields['entropy'][..., k]; h = fields['specific_internal'][..., k]
        ha, sa = float(h.flat[0]), float(s.flat[0]); dh, ds = h-ha, s-sa
        phi = operator.potential.apply(gamma*h-.5*s*h-.5*ha*s+.5*sa*h)
        remainder = -.5*(operator.mapping.scalar(dh)[..., None]*operator.gradient(ds)
                         - operator.mapping.scalar(ds)[..., None]*operator.gradient(dh))
        reduced = a/E*weak_vector(geo, remainder); potential = a/E*geo.gradient(phi)
        full = m.mass@potential+reduced
        if fingerprint(full) != fingerprint(original['species_forcing'][k]):
            raise ValueError('Original per-species force must match bitwise')
        result.append((reduced, potential, full))
    for index, key in ((0, 'reduced_forcing'), (1, 'potential_electric')):
        if fingerprint(sum(item[index] for item in result)) != fingerprint(original[key]):
            raise ValueError('Original complete species split must match bitwise')
    return result


def decompose(system, state, field, original):
    d, geo = system.fluid, system.geometry
    model = RegularCurrent(system.mapping, state[..., 0], system.ion_mass)
    bulk, _ = d.derivative(state)
    unconstrained, _ = RegularCoordinateFluid.derivative(d, state)
    bulk_reaction = np.zeros_like(state); bulk_reaction[-1, -1, :, 1] = -unconstrained[-1, -1, :, 1]
    if not np.array_equal(unconstrained+bulk_reaction, bulk):
        raise ValueError('Explicit original bulk wall reaction required')
    em = model.electromagnetic_pair(state, field, bulk[..., 0])
    transport_op, species_op, induction = system.operators(model)
    transport = transport_op.pair(state, em['current'], bulk[..., 0])
    species = species_op.pair(state, em['current'])
    separated = species_parts(species_op, state, species)
    assembled = bulk+em['fluid_rate']+transport['fluid_rate']+species['fluid_rate']
    reaction = np.zeros_like(state); reaction[-1, -1, :, 1] = -assembled[-1, -1, :, 1]
    assembled[-1, -1, :, 1] = 0.
    if fingerprint(assembled) != fingerprint(original.fluid): raise ValueError('Original fluid assembly required')
    forcing = em['forcing']+transport['forcing']+species['reduced_forcing']
    full_forcing = em['forcing']+transport['forcing']+species['forcing']
    zero_e = np.zeros(model.mass.shape[0]); zero_b = np.zeros(geo.hessian.shape[0])
    arrays = dict(bulk_before_wall=unconstrained, bulk_wall_reaction=bulk_reaction,
                  original_reduced_forcing=forcing, original_full_forcing=full_forcing,
                  original_current=em['current'])
    fluid = {name: np.zeros_like(state) for name in GROUPS}
    fluid['bulk'] = bulk; fluid['ideal_lorentz'] = em['fluid_rate']
    fluid['stretching_drift_stress'] = transport['fluid_rate']; fluid['wall_reaction'] = reaction
    fluid['electron_species'][..., 4] = species['fluid_rate'][..., 4]
    fluid['ion_species'][..., 5] = species['fluid_rate'][..., 5]
    forces = dict(ideal_lorentz=em['ideal_forcing'], hall=em['hall_forcing'],
                  density_inertia=em['density_forcing'], current_advection=transport['advection'],
                  stretching_drift_stress=transport['stretching'], self_drift=transport['self_drift'],
                  electron_species=separated[0][0], ion_species=separated[1][0])
    potential = dict(electron_species=separated[0][1], ion_species=separated[1][1])
    records = []
    for name in GROUPS:
        f = forces.get(name, zero_e); ep = potential.get(name, zero_e)
        full = model.mass@ep+f
        arrays.update({name+'_'+key: value for key, value in dict(fluid=fluid[name], forcing=f,
                       full_forcing=full, potential_electric=ep).items()})
        row = dict(name=name, status='running'); records.append(row)
        if name not in forces:
            b, j, electric = zero_b.copy(), zero_e.copy(), zero_e.copy()
            row.update(status='complete', residual=0., original_full_ohm_backward=0., residual_components={})
        else:
            try: solved = induction.solve(f)
            except (ValueError, RuntimeError) as error:
                row.update(status='induction_rejected', exception_type=type(error).__name__,
                           message=str(error), attempts=getattr(induction, 'attempt_trace', []))
                continue
            b, j, electric = solved.magnetic, solved.current, solved.electric+ep
            lhs = model.mass@electric; inertia = model.inertia@j
            restored = backward(lhs-inertia-full, abs(lhs)+abs(inertia)+abs(full))
            if restored is None or max(solved.residual, restored) > system.solver_options['tolerance']:
                row.update(status='restored_ohm_rejected', residual=solved.residual,
                           original_full_ohm_backward=restored, residual_components=solved.residual_components)
                continue
            row.update(status='complete', residual=solved.residual, original_full_ohm_backward=restored,
                       iterations=solved.iterations, residual_components=solved.residual_components)
        arrays.update({name+'_'+key: value for key, value in dict(magnetic=b, current=j, electric=electric).items()})
    result = dict(groups=records, all_parts_complete=all(r['status']=='complete' for r in records),
                  forcing_sum=difference(sum(forces.values()), forcing),
                  full_forcing_sum=difference(sum(arrays[n+'_full_forcing'] for n in GROUPS), full_forcing))
    if result['all_parts_complete']:
        total = {key: sum(arrays[n+'_'+key] for n in GROUPS) for key in ('fluid', 'magnetic', 'current', 'electric')}
        result['sum_vs_original'] = {key: difference(value, getattr(original, key)) for key, value in total.items()}
        lhs = model.mass@total['electric']; inert = model.inertia@total['current']
        amp = model.mass@total['current']; weak_amp = -geo.faraday.T@(geo.hessian@total['magnetic'])
        curl = geo.project_curl(difference_curl(geo, total['electric']))
        checks = dict(ohm=backward(lhs-inert-full_forcing, abs(lhs)+abs(inert)+abs(full_forcing)),
                      ampere=backward(amp-weak_amp, abs(amp)+abs(weak_amp)),
                      faraday=backward(curl-total['magnetic'], abs(geo.faraday)@abs(total['electric'])+abs(total['magnetic'])))
        result.update(recombined_original_backward=checks,
                      recombined_original_bounds_pass=all(v is not None and v<=system.solver_options['tolerance'] for v in checks.values()))
        arrays.update({'sum_'+key: value for key, value in total.items()})
    return result, arrays


def run(destination):
    destination = Path(destination); probe = Path('outputs/trajectory_tangent_probe')
    if destination.exists(): raise ValueError('Fresh decomposition destination required')
    verification = Path('outputs/trajectory_tangent_probe_replay_verification/summary.json')
    replay = json.loads(verification.read_text()); original = json.loads((probe/'summary.json').read_text())
    prior_source = json.loads((probe/'provenance.json').read_text())
    if not original['terminal'] or not original['all_cases_complete'] or not replay['terminal'] or not replay['bitwise_metrics']:
        raise ValueError('Complete original tangent replay required')
    inputs = {str(verification): sha(verification), str(probe/'provenance.json'): sha(probe/'provenance.json')}
    for side in ('original', 'repeated'):
        root = Path(replay[side+'_directory']); p = root/'summary.json'
        if sha(p) != replay[side+'_summary_sha256']: raise ValueError('Original replay summaries required')
        inputs[str(p)] = sha(p)
        for name, row in replay['arrays'].items():
            if row['original_sha256'] != row['repeated_sha256'] or sha(root/name) != row[side+'_sha256']:
                raise ValueError('Every original raw replay file required')
            inputs[str(root/name)] = sha(root/name)
    if any(sha(p)!=h for p,h in prior_source['source_files'].items()): raise ValueError('Original probe source required')
    identity = provenance(['examples.trajectory_spatial_parts']); destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=False, source_sha256=identity['source_sha256'], input_sha256=inputs,
                  probe=str(probe), spec=original['spec'], groups=list(GROUPS), cases=[])
    def save():
        p = destination/'summary.pending.json'; p.write_text(json.dumps(result, indent=2)+'\n'); p.replace(destination/'summary.json')
    save()
    try:
        _, _, base, system, bias, _, _ = fixture(original['spec'])
        if system.solver_options['tolerance'] != 1e-10: raise ValueError('Original bound required')
        coordinates = ModalWallCoordinates(system, base, bias)
        with np.load(probe/'context.npz') as raw: context = dict(raw)
        if fingerprint(base) != fingerprint(context['state']): raise ValueError('Original base required')
        base_rate = system.derivative(base, bias)
        with np.load(probe/'base_rate.npz') as raw:
            if any(fingerprint(getattr(base_rate, k))!=fingerprint(raw[k+'_rate' if k!='electric' else k])
                   for k in ('fluid','magnetic','current','electric')): raise ValueError('Original base RHS required')
        for index, case in enumerate(original['cases']):
            row = dict(index=index, direction=case['direction'], amplitude=case['amplitude'], signs=[])
            result['cases'].append(row); save(); all_arrays={}
            with np.load(probe/case['input_file']) as raw: states = dict(raw)
            for sign_row in case['signs']:
                sign = sign_row['sign']; state = states[sign+'_state']; field = coordinates.field(states[sign+'_magnetic'])
                before = (fingerprint(state), fingerprint(field.packed_increment(system.geometry)))
                full = system.derivative(state, field)
                with np.load(probe/sign_row['arrays_file']) as raw:
                    matches = {k: fingerprint(getattr(full, k))==fingerprint(raw[k+'_rate' if k!='electric' else k])
                               for k in ('fluid','magnetic','current','electric')}
                if not all(matches.values()): raise ValueError('Complete signed original RHS must match bitwise')
                report, arrays = decompose(system, state, field, full)
                unchanged = before==(fingerprint(state), fingerprint(field.packed_increment(system.geometry)))
                if not unchanged: raise ValueError('Actual decomposition input modified')
                arrays.update(state=state, magnetic_input=states[sign+'_magnetic'])
                path = destination/f'case_{index:02}_{sign}.npz'; np.savez_compressed(path, **arrays)
                report.update(sign=sign, arrays_file=path.name, arrays_sha256=sha(path),
                              original_full_rhs_bitwise=matches, original_residual=full.residual,
                              actual_inputs_byte_unchanged=unchanged)
                row['signs'].append(report); all_arrays[sign]=arrays; save()
                print('spatial parts', index, sign, report['all_parts_complete'], report.get('recombined_original_bounds_pass'), flush=True)
            row['all_parts_complete'] = all(s['all_parts_complete'] for s in row['signs'])
            if row['all_parts_complete']:
                direction = context[case['direction']+'_coordinate']; denominator = float(direction@direction)
                projections={}
                for name in GROUPS:
                    derivatives={key:(all_arrays['plus'][name+'_'+key]-all_arrays['minus'][name+'_'+key])/(2*case['amplitude'])
                                 for key in ('fluid','magnetic')}
                    action=coordinates.pack(**derivatives)
                    projections[name]=float(direction@action/denominator)
                row['coordinate_projection_contributions']=projections
                row['original_coordinate_projection']=case['coordinate_projection']
                row['coordinate_projection_sum_difference']=sum(projections.values())-case['coordinate_projection']
            save()
        result.update(all_cases_complete=all(c['all_parts_complete'] for c in result['cases']),
                      all_recombined_original_bounds_pass=all(s.get('recombined_original_bounds_pass',False) for c in result['cases'] for s in c['signs']))
    except BaseException as error:
        result.update(exception_type=type(error).__name__,exception_message=str(error)); raise
    finally:
        result['terminal']=True
        if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.trajectory_spatial_parts'])!=identity:
            result['source_or_input_changed']=True; save(); raise ValueError('Source or input changed')
        save()
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
