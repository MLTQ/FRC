"""Isolate fast saved response and test the actual stationary spatial RHS."""
import json
from pathlib import Path
import numpy as np
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.uniform_wall_plasma import preconditioned
from .matched_mode_comparison import load
from .matched_mode_fixture import fixture
from .current_transport_sampling import TransportSampler
from .extracted_response_fields import fields, inner, localization
from .matched_mode_rates import fingerprint
from .implicit_plasma_records import provenance, freeze, sha


def calculate(audit_directory):
    audit, _ = load(audit_directory)
    directories = [Path(p) for p in audit['input_directories']]
    records = [json.loads((p/'summary.json').read_text()) for p in directories]
    mode, d, base, original, bias, _, _ = fixture(audit['spec'])
    system = preconditioned(original); c = ModalWallCoordinates(system, base, bias)
    if audit['steps'] % 5: raise ValueError('Five equal time segments required')
    snapshots = {}; hashes = {}
    for index in (2*audit['steps']//5, 3*audit['steps']//5, 4*audit['steps']//5, audit['steps']):
        states = []
        for directory, record in zip(directories, records):
            entry = record['records'][index-1]; raw = directory/entry['arrays_file']
            if sha(raw) != entry['arrays_sha256']: raise ValueError('Unchanged checkpoint required')
            hashes[str(raw)] = sha(raw); states.append(dict(np.load(raw)))
        snapshots[index] = c.pack(states[1]['state']-states[0]['state'],
            states[1]['magnetic']-states[0]['magnetic'])/(2*audit['epsilon'])
    segment = audit['steps']//5; interval = audit['duration']/5
    slow_growth = np.exp(mode.growth*2*interval)
    early = snapshots[4*segment]-slow_growth*snapshots[2*segment]
    late = snapshots[5*segment]-slow_growth*snapshots[3*segment]
    scale = float(np.max(abs(late))); direction = late/scale
    du, db = c.unpack(direction); rates = []; diagnostics = []
    for epsilon in (1e-4, 1e-5):
        pair = [system.derivative(base+sign*epsilon*du, c.field(sign*epsilon*db)) for sign in (-1, 1)]
        value = c.pack(pair[1].fluid-pair[0].fluid, pair[1].magnetic-pair[0].magnetic)/(2*epsilon)
        rates.append(value)
        diagnostics.append(dict(epsilon=epsilon, residuals=[p.residual for p in pair],
            induction=[p.diagnostics['induction_residuals'] for p in pair]))
    arrays = dict(base=base, early_direction=early, late_direction=late, normalized_direction=direction,
        fluid_direction=du, magnetic_direction=db, rate_eps1e4=rates[0], rate_eps1e5=rates[1])
    quadrature = []
    for order in (16, 20):
        sample = TransportSampler(d, system.geometry.mesh, order)
        previous = fields(sample, system, base, *c.unpack(early))
        current = fields(sample, system, base, *c.unpack(late))
        tangent = [fields(sample, system, base, *c.unpack(rate*scale)) for rate in rates]
        channels = {}; prints = {}
        for name, value in current.items():
            a = previous[name]; norm = inner(sample, value, value)
            gain = inner(sample, a, value)/inner(sample, a, a)
            measured = np.log(gain)/interval if gain > 0 else None
            result = dict(time_window_gain=gain, inferred_growth=measured,
                window_shape_error=float(np.sqrt(inner(sample, value-gain*a, value-gain*a)/norm)))
            loc, raw = localization(sample, value); result.update(loc)
            for key, item in raw.items(): arrays[f'q{order}_{name}_{key}'] = item
            result['rhs_controls'] = []
            for rate in tangent:
                rayleigh = inner(sample, value, rate[name])/norm
                residual = rate[name]-rayleigh*value
                result['rhs_controls'].append(dict(projected_growth=rayleigh,
                    off_direction_rate_relative=float(np.sqrt(inner(sample, residual, residual)/inner(sample, rate[name], rate[name]))),
                    inferred_growth_residual_relative=None if measured is None else float(np.sqrt(
                        inner(sample, rate[name]-measured*value, rate[name]-measured*value)/inner(sample, rate[name], rate[name])))))
            diff = tangent[1][name]-tangent[0][name]
            result['amplitude_rhs_relative_difference'] = float(np.sqrt(inner(sample, diff, diff)/inner(sample, tangent[1][name], tangent[1][name])))
            channels[name] = result; prints[name] = fingerprint(value)
        arrays[f'q{order}_radius'] = sample.r
        quadrature.append(dict(order=order, channels=channels, physical_fingerprints=prints))
    return dict(audit_directory=str(audit_directory), input_summary_sha256=sha(Path(audit_directory)/'summary.json'),
        input_checkpoints=hashes, spec=audit['spec'], duration=audit['duration'], epsilon=audit['epsilon'],
        stationary_continuum_growth=mode.growth, direction_scale=scale, rhs_diagnostics=diagnostics,
        quadrature=quadrature), arrays


def run(audit_directory, destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    identity = provenance(['examples.extracted_response_probe']); destination.mkdir(parents=True); freeze(destination, identity)
    result, arrays = calculate(audit_directory)
    raw = destination/'response.npz'; np.savez_compressed(raw, **arrays)
    if provenance(['examples.extracted_response_probe']) != identity: raise RuntimeError('Source changed')
    result.update(terminal=True, source_sha256=identity['source_sha256'], arrays_file=raw.name, arrays_sha256=sha(raw))
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print('complete extracted response', destination, flush=True)
    return result
