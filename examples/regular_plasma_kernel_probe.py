"""Capture the original strict cancellation failure without accepting a bad rate."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from plasma_bench.kernel_induction import KernelInduction, separate, difference_curl
from plasma_bench.inertial_induction import InertialInduction
from plasma_bench.regular_plasma import RegularPlasma
from .regular_plasma_spectrum import run_case
from .regular_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing to overwrite failure capture')
    identity = provenance()
    root = Path(__file__).resolve().parents[1]
    for name in ('examples/regular_plasma_kernel_probe.py', 'examples/regular_plasma_spectrum.py'):
        identity['source_files'][name] = sha(root/name)
    identity['source_sha256'] = hashlib.sha256(json.dumps(identity['source_files'], sort_keys=True).encode()).hexdigest()
    destination.mkdir(parents=True); freeze(destination, identity)
    original = KernelInduction.solve
    original_derivative = RegularPlasma.derivative
    captured = {}
    inputs = {}

    def observe_input(self, state, field):
        inputs.update(fluid_state=np.array(state, copy=True),
                      magnetic_state=field.packed_increment(self.geometry).copy())
        return original_derivative(self, state, field)

    def observed(self, forcing):
        try:
            return original(self, forcing)
        except RuntimeError as exc:
            electric0 = self.mass_inverse(forcing)
            small, kernel = separate(self.model.geometry, electric0)
            result = InertialInduction.solve(self, self.model.mass@small)
            full = result.electric+kernel
            curl = self.model.geometry.faraday
            captured.update(inputs, forcing=forcing, electric0=electric0, small=small, kernel=kernel,
                solved_electric=result.electric, restored_electric=full, magnetic_rate=result.magnetic,
                matrix_curl=curl@full, difference_curl=difference_curl(self.model.geometry, full),
                roundoff_bound=abs(curl)@abs(full)+abs(result.magnetic),
                current_rate=result.current, number_density=self.model.density_map.point_values)
            np.savez_compressed(destination/'failure.npz', **captured)
            norm = np.linalg.norm
            report = dict(source_sha256=identity['source_sha256'], status='rejected_as_required', error=repr(exc),
                arrays_sha256=sha(destination/'failure.npz'),
                norms={k: float(norm(v)) for k, v in captured.items()},
                difference_vs_matrix=float(norm(captured['difference_curl']-captured['matrix_curl'])),
                matrix_vs_rate=float(norm(captured['matrix_curl']-result.magnetic)),
                difference_vs_rate=float(norm(captured['difference_curl']-result.magnetic)),
                inner_residual=result.residual, inner_residuals=result.residual_components)
            report['component_spreads'] = np.ptp(inputs['fluid_state'], axis=(0, 1, 2)).tolist()
            (destination/'failure.json').write_text(json.dumps(report, indent=2)+'\n')
            raise

    KernelInduction.solve = observed
    RegularPlasma.derivative = observe_input
    try:
        run_case(dict(elements=1, points=8, field_nr=6, field_nz=8, epsilon=1e-4))
    except RuntimeError:
        if not captured:
            raise
    finally:
        KernelInduction.solve = original
        RegularPlasma.derivative = original_derivative
    if not captured:
        raise RuntimeError('Expected original rejection was not reproduced')
    assert all(sha(root/name) == value for name, value in identity['source_files'].items())
    print('Original rejection captured; no failed rate accepted', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    run(parser.parse_args().destination)
