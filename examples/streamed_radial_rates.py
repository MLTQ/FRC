"""Full signed spatial workflow with bounded-memory audits and early records."""
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.regular_plasma_step import field_from_packed
from plasma_bench.split_inverse_wall_plasma import split_inverse
from .archived_contact_mode import ArchivedContactMode
from .radial_contact_fixture import fixture
from .streamed_rate_audit import measure
from .streamed_interface_audit import measure as interfaces
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, specs, reference_directory='outputs/radial16_reference', reference_index=1, chunk_cells=2):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing overwrite')
    if type(chunk_cells) is not int or chunk_cells < 1:
        raise ValueError('Positive integer radial cell batch required')
    mode = ArchivedContactMode(reference_directory, reference_index)
    identity = provenance(['examples.streamed_radial_rates'])
    destination.mkdir(parents=True); freeze(destination, identity)
    records = []
    summary = dict(terminal=False, source_sha256=identity['source_sha256'], cases=records,
        reference_directory=str(reference_directory), reference_index=reference_index, chunk_cells=chunk_cells)
    def save():
        (destination/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    save()
    for spec in specs:
        start = time.monotonic(); epsilon = spec.get('epsilon', 1e-4)
        if not np.isfinite(epsilon) or epsilon <= 0:
            raise ValueError('Positive finite perturbation amplitude required')
        index = len(records)
        record = dict(spec=spec, epsilon=epsilon, input_sha256=mode.input_sha256,
            input_summary_sha256=mode.input_summary_sha256, radial_mode=mode.radial_mode,
            reference_points=mode.points, continuum_growth=mode.growth, status='initializing')
        records.append(record); save()
        print('streamed case', index, spec, 'initializing', flush=True)
        d, u, system, bias, du, db = fixture(spec, mode); system = split_inverse(system)
        arrays = dict(state=u, fluid_perturbation=du, magnetic_perturbation=db); rates = []
        path = destination/f'case_{index:02}.npz'
        def save_arrays():
            np.savez_compressed(path, **arrays)
            record.update(arrays_file=path.name, arrays_sha256=sha(path))
            save()
        for name, sign in [('minus', -1), ('plus', 1)]:
            state = u+sign*epsilon*du
            field = field_from_packed(system.geometry, sign*epsilon*db, bias.axial_bias)
            arrays.update({name+'_state': state, name+'_magnetic': field.packed_increment(system.geometry)})
            record.update(status='evaluating_'+name); save_arrays()
            print('streamed case', index, name, 'full RHS', flush=True)
            try:
                rate = system.derivative(state, field)
            except (RuntimeError, ValueError) as error:
                record.update(status='failed', failed_sign=name, failure_type=type(error).__name__, failure_message=str(error))
                break
            rates.append(rate)
            arrays.update({name+'_'+key: value for key, value in dict(fluid_rate=rate.fluid,
                magnetic_rate=rate.magnetic, electric=rate.electric, current_rate=rate.current).items()})
            save_arrays()
        if len(rates) == 2:
            record.update(status='auditing', quadrature=[], max_induction_residual=max(r.residual for r in rates),
                rate_diagnostics=[r.diagnostics for r in rates]); save()
            for order in (16, 20):
                print('streamed case', index, 'physical quadrature', order, flush=True)
                physical = measure(mode, d, system, u, du, db, epsilon, rates, order, chunk_cells)
                face, raw = interfaces(mode, d, system, u, db, order, chunk_cells)
                physical['representation_and_faces'] = face
                arrays.update({f'q{order}_{key}': value for key, value in raw.items()})
                record['quadrature'].append(physical); save_arrays()
            record['status'] = 'complete'
        record['elapsed_seconds'] = time.monotonic()-start
        (destination/f'case_{index:02}.json').write_text(json.dumps(record, indent=2)+'\n'); save()
        print('streamed case', index, record['status'], 'seconds', record['elapsed_seconds'], flush=True)
    if provenance(['examples.streamed_radial_rates']) != identity:
        raise RuntimeError('Source changed')
    summary.update(terminal=True, all_cases_complete=all(r['status'] == 'complete' for r in records)); save()
    return summary
