"""Observe actual assembled operators during complete original/tensor pilots."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.sparse import issparse
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.tensor_mass_geometry import TensorMassGeometry
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from .direct_metric_rates import run as original_rates
from .tensor_metric_rates import run as tensor_rates
from .wall_diagnostic_replay import verify
from .implicit_plasma_records import provenance, freeze, sha


def fingerprint(value):
    if issparse(value):
        return dict(format=value.format, shape=list(value.shape),
                    data=fingerprint(value.data), indices=fingerprint(value.indices),
                    indptr=fingerprint(value.indptr))
    value = np.asarray(value)
    return dict(shape=list(value.shape), dtype=str(value.dtype),
                sha256=hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest())


def operators(geometry):
    names = ('radial_points', 'axial_points', 'wr', 'wz', 'Nx', 'dNx', 'Ex', 'Dx',
             'Nz', 'dNz', 'Ez', 'Dz', 'radial_e', 'theta_e', 'radial_b', 'offsets',
             'hessian', 'faraday', 'r', 'z', 'weight', 'wall_trace',
             'radial_projection', 'projection')
    result = {name: fingerprint(getattr(geometry, name)) for name in names if hasattr(geometry, name)}
    result.update({f'mass_{k}': fingerprint(m) for k, m in enumerate(geometry.mass)})
    return result


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Fresh observation directory required')
    identity = provenance(['examples.tensor_operator_workflow', 'examples.streamed_radial_rates',
                           'examples.cached_streamed_rates', 'examples.streamed_interface_audit'])
    destination.mkdir(parents=True); freeze(destination, identity)
    records = {}; comparisons = {}
    spec = dict(elements=4, points=64, field_nr=16, field_nz=64, epsilon=1e-4, tolerance=1e-10)
    for label, cls, calculate, archived in (
            ('original', CompatibleGeometry, original_rates, 'outputs/tensor_mass_pilot_original'),
            ('tensor', TensorMassGeometry, tensor_rates, 'outputs/tensor_mass_pilot')):
        captured = []; originals = (cls.__init__, StrongWallGeometry.__init__)
        def observe(init, stage):
            def initialize(self, *args, **kwargs):
                init(self, *args, **kwargs)
                captured.append(dict(stage=stage, operators=operators(self)))
            return initialize
        cls.__init__ = observe(originals[0], 'base')
        StrongWallGeometry.__init__ = observe(originals[1], 'strong_wall')
        target = destination/label
        try:
            result = calculate(target, [spec], reference_directory='outputs/stationary_contact_audit', reference_index=3)
        finally:
            cls.__init__, StrongWallGeometry.__init__ = originals
        if not result['terminal'] or not result['all_cases_complete']:
            raise ValueError('Complete signed physical workflow required')
        if [row['stage'] for row in captured] != ['base', 'strong_wall']:
            raise ValueError('Both actual construction stages required')
        record = destination/(label+'_operators.json')
        record.write_text(json.dumps(captured, indent=2)+'\n'); records[label] = captured
        verify(archived, target, destination/(label+'_verification'))
        comparisons[label] = dict(operators_file=record.name, operators_sha256=sha(record),
                                 full_archived_workflow_replayed=True)
    if records['original'] != records['tensor']:
        raise ValueError('Actual assembled operators changed')
    if provenance(['examples.tensor_operator_workflow', 'examples.streamed_radial_rates',
                   'examples.cached_streamed_rates', 'examples.streamed_interface_audit']) != identity:
        raise ValueError('Source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], spec=spec,
                  operator_fingerprints_identical=True, cases=comparisons,
                  scope='Read-only observations during complete signed coarse physical workflows. Both independently replay archived inputs, rates and audits; matrix identity here is measured on this coarse grid only. Inverse results are not claimed bitwise equal between algorithms.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
