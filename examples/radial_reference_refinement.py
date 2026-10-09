"""Independent axial convergence at one fixed radial Bessel wavenumber."""
import json
from pathlib import Path
import time
import numpy as np
from .stationary_contact_continuum import operator
from .implicit_plasma_records import provenance,freeze,sha


def run(destination,points=(64,128,256),radial_mode=16):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    identity=provenance(['examples.radial_reference_refinement']);destination.mkdir(parents=True);freeze(destination,identity)
    continuum=[]
    for count in points:
        start=time.monotonic();matrix,arrays=operator(count,radial_mode)
        values,vectors=np.linalg.eig(matrix);index=int(values.real.argmax());vector=vectors[:,index]
        physical=(vector*arrays['scale']).reshape(9,count)
        divergence=arrays['q']*physical[6]+arrays['derivative']@physical[8]
        denominator=np.linalg.norm(arrays['q']*physical[6])+np.linalg.norm(arrays['derivative']@physical[8])
        path=destination/f'continuum_{len(continuum):02}.npz'
        np.savez_compressed(path,**arrays,operator=matrix,eigenvalues=values,leading_vector=vector)
        record=dict(points=count,radial_mode=radial_mode,q=float(arrays['q']),maximum_real=float(values.real[index]),
            leading_imaginary=float(values.imag[index]),divergence_relative=float(np.linalg.norm(divergence)/denominator),
            eigenpair_backward=float(np.linalg.norm(matrix@vector-values[index]*vector)/(np.linalg.norm(matrix)*np.linalg.norm(vector))),
            arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start)
        continuum.append(record)
        print('independent radial reference',radial_mode,count,record['maximum_real'],flush=True)
    if provenance(['examples.radial_reference_refinement'])!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],continuum=continuum)
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
