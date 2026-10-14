"""Continuum prediction of the added numerical term on the resolved contact."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.linalg import expm
from .stationary_contact_continuum import operator
from .continuum_contact_mode import ContinuumContactMode
from .implicit_plasma_records import provenance, freeze


def run(destination, eta6, duration=1e-4):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Fresh reference output required')
    if not np.isfinite([eta6, duration]).all() or eta6 < 0 or duration <= 0:
        raise ValueError('Nonnegative eta6 and positive duration required')
    mode = ContinuumContactMode()
    matrix, arrays = operator(mode.points)
    curl = arrays['curl']
    c2 = curl@curl
    damping = -eta6/MU0*np.linalg.solve(arrays['inertia'], c2@c2@c2)
    magnetic_scale = arrays['scale'][6*mode.points:]
    modified = matrix.copy()
    modified[6*mode.points:, 6*mode.points:] += damping*magnetic_scale[None, :]/magnetic_scale[:, None]
    initial = mode.physical_nodes.ravel()/arrays['scale']
    original = expm(duration*matrix)@initial
    evolved = expm(duration*modified)@initial
    target = initial*np.exp(mode.growth*duration)
    norms = lambda v: np.linalg.norm(v.reshape(9, mode.points), axis=1)
    identity = provenance(['examples.hyperdiffusive_contact_reference'])
    destination.mkdir(parents=True); freeze(destination, identity)
    np.savez_compressed(destination/'reference.npz', initial=initial, original=original,
        modified=evolved, target=target, original_operator=matrix, modified_operator=modified, **arrays)
    record = dict(terminal=True, eta6_ohm_m5=eta6, duration=duration,
        source_sha256=identity['source_sha256'], input_sha256=mode.input_sha256,
        names=['density','vr','vtheta','vz','electron_entropy','ion_entropy','Br','Btheta','Bz'],
        numerical_change_relative_by_component=(norms(evolved-original)/norms(original)).tolist(),
        original_exponential_error_by_component=(norms(original-target)/norms(target)).tolist(),
        scope='Independent Fourier/Bessel linear continuum evolution, including prescribed eta6 curl^4 J. Numerical heat is quadratic and absent at linear order. Amplitude norms share each component radial Bessel factor; no nonlinear or discrete spatial trajectory is used.')
    (destination/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2), flush=True)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    parser.add_argument('--eta6', required=True, type=float)
    args = parser.parse_args(); run(args.destination, args.eta6)
