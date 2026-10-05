"""Frozen full-coordinate and trajectory checks of the opt-in refined assembly."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky, solve_triangular
import examples.regular_plasma_reference as fixtures
from .regular_plasma_spectrum import run_case as coordinates
from .regular_plasma_evolution import run_case as evolve
from .regular_plasma_records import provenance, freeze, sha
from plasma_bench.refined_regular_plasma import RefinedRegularPlasma
from plasma_bench.refined_kernel_induction import RefinedKernelInduction


def run(source, destination):
    source = Path(source); destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing to overwrite refinement study')
    identity = provenance(); root = Path(__file__).resolve().parents[1]
    assert json.loads((source/'provenance.json').read_text()) == identity
    for name in ('examples/regular_plasma_refinement_study.py', 'examples/regular_plasma_spectrum.py',
                 'plasma_bench/refined_kernel_induction.py', 'plasma_bench/refined_regular_plasma.py',
                 'tests/test_refined_kernel.py'):
        identity['source_files'][name] = sha(root/name)
    identity['source_sha256'] = hashlib.sha256(json.dumps(identity['source_files'], sort_keys=True).encode()).hexdigest()
    destination.mkdir(parents=True); freeze(destination, identity)
    original_class = fixtures.RegularPlasma; original_solve = RefinedKernelInduction.solve
    stats = []

    def observed(self, forcing):
        result = original_solve(self, forcing)
        stats.append((result.residual, result.residual_components['outer_refinements']))
        return result

    fixtures.RegularPlasma = RefinedRegularPlasma; RefinedKernelInduction.solve = observed
    try:
        target_record = json.loads((source/'case_20.json').read_text())
        path = source/target_record['arrays_file']; assert sha(path) == target_record['arrays_sha256']
        target = dict(np.load(path)); root_metric = cholesky(target['energy_metric'], lower=False)
        reference_norm = np.linalg.norm(target['energy_weighted_operator'])
        rows = []
        for i, epsilon in enumerate((1e-4, 1e-5)):
            spec = dict(elements=1, points=8, field_nr=6, field_nz=8, epsilon=epsilon)
            stats.clear(); result, arrays = coordinates(spec)
            result['solve_count'] = len(stats)
            result['refined_solve_count'] = sum(n>0 for _, n in stats)
            result['max_outer_refinements'] = max(n for _, n in stats)
            result['max_original_residual'] = max(r for r, _ in stats)
            difference = root_metric@(arrays['operator']-target['operator'])
            difference = solve_triangular(root_metric.T, difference.T, lower=True).T
            result['independent_linear_operator_relative_error'] = float(np.linalg.norm(difference)/reference_norm)
            repeated, replay = coordinates(spec)
            assert all(np.array_equal(arrays[k], replay[k]) for k in arrays)
            for key, value in repeated.items():
                if key != 'elapsed_seconds': assert result[key] == value
            result['bitwise_replayed'] = True
            raw = destination/f'case_{i:02d}.npz'; np.savez_compressed(raw, **arrays)
            result.update(source_sha256=identity['source_sha256'], arrays_sha256=sha(raw), arrays_file=raw.name)
            raw.with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n'); rows.append(result)
            print(raw.name, 'all coordinate perturbations verified', flush=True)
        spec = dict(kind='wave', elements=4, amplitude=1., steps=120)
        stats.clear(); result, arrays = evolve(spec)
        old = json.loads((source/'case_12.json').read_text()); path = source/old['arrays_file']
        assert sha(path) == old['arrays_sha256']; target_arrays = dict(np.load(path))
        assert all(np.array_equal(arrays[k], target_arrays[k]) for k in arrays)
        assert all(old[k] == result[k] for k in result if k != 'elapsed_seconds')
        result.update(legacy_full_trajectory_bitwise=True, refined_solve_count=sum(n>0 for _, n in stats),
                      max_original_residual=max(r for r, _ in stats))
        raw = destination/'case_02.npz'; np.savez_compressed(raw, **arrays)
        result.update(source_sha256=identity['source_sha256'], arrays_sha256=sha(raw), arrays_file=raw.name)
        raw.with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n'); rows.append(result)
        print(raw.name, 'complete wave trajectory identical to original', flush=True)
    finally:
        fixtures.RegularPlasma = original_class; RefinedKernelInduction.solve = original_solve
    assert all(sha(root/name) == value for name, value in identity['source_files'].items())
    (destination/'summary.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'], terminal=True,
                                                          records=rows), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('source'); parser.add_argument('destination')
    args = parser.parse_args(); run(args.source, args.destination)
