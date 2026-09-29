"""Volume-weighted field comparisons for the recorded driven boundary study."""
from pathlib import Path
import hashlib
import json
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh


def read_case(directory, name):
    directory = Path(directory)
    row = json.loads((directory/f"{name}.json").read_text())
    if not row["completed"]:
        raise ValueError(f"{name} has not completed")
    path = directory/row["fields"]["file"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != row["fields"]["sha256"]:
        raise ValueError(f"{name} fields do not match their recorded hash")
    with np.load(path, allow_pickle=False) as archive:
        fields = {key: archive[key] for key in archive.files}
    return row, fields


def restrict_cells(values, fine, coarse):
    """Exact cylindrical volume average from nested uniform fine cells."""
    factor = fine.nr//coarse.nr
    if (factor < 1 or fine.nr != factor*coarse.nr or fine.nz != factor*coarse.nz or
        fine.radius != coarse.radius or fine.half_length != coarse.half_length):
        raise ValueError("Field comparisons require nested meshes of equal extent")
    tail = values.shape[2:]
    weight = fine.volumes.reshape((*fine.volumes.shape,)+(1,)*len(tail))
    blocks = (values*weight).reshape(coarse.nr, factor, coarse.nz, factor, *tail)
    result = blocks.sum(axis=(1, 3))/coarse.volumes.reshape((*coarse.volumes.shape,)+(1,)*len(tail))
    np.testing.assert_allclose(np.sum(result*coarse.volumes.reshape((*coarse.volumes.shape,)+(1,)*len(tail)), axis=(0, 1)),
                               np.sum(values*weight, axis=(0, 1)), rtol=2e-12, atol=1e-14)
    return result


def compare(a, b):
    """Differences normalized by the finer reference's physical perturbations."""
    ar, af = a; br, bf = b
    am, bm = MHDMesh(**ar["fluid_mesh"]), MHDMesh(**br["fluid_mesh"])
    factor = bm.nr//am.nr
    rows = []
    for item in ar["fields"]["times"]:
        prefix = item["prefix"]
        q = af[prefix+"_primitive"]
        ref = restrict_cells(bf[prefix+"_primitive"], bm, am)
        def error(value, exact):
            dims = (1,)*max(value.ndim-2, 0)
            weight = am.volumes.reshape(am.volumes.shape+dims)
            return float(np.sqrt(np.sum((value-exact)**2*weight)/np.sum(exact**2*weight)))
        bias = ar["source_configuration"]["bias_t"]
        current = af[prefix+"_current"]
        jr = bf[prefix+"_current"][::factor, ::factor]
        node_weight = am.rf[:, None]*np.ones((am.nr+1, am.nz+1))
        node_weight[-1] *= .5; node_weight[:, [0, -1]] *= .5
        rows.append({"time_s": item["time_s"],
                     "velocity_relative_L2": error(q[..., 1:3], ref[..., 1:3]),
                     "density_perturbation_relative_L2": error(q[..., 0]-1, ref[..., 0]-1),
                     "pressure_perturbation_relative_L2": error(q[..., 3]-1, ref[..., 3]-1),
                     "magnetic_perturbation_relative_L2": error(q[..., 4:6]-[0, bias], ref[..., 4:6]-[0, bias]),
                     "current_relative_L2": float(np.sqrt(np.sum((current-jr)**2*node_weight)/np.sum(jr**2*node_weight)))})
    return {"coarse_or_test": ar["name"], "reference": br["name"], "history": rows}


def summarize(case):
    row, fields = case
    mesh = MHDMesh(**row["fluid_mesh"])
    samples = []
    for item in row["fields"]["times"]:
        prefix = item["prefix"]
        q, j = fields[prefix+"_primitive"], fields[prefix+"_current"]
        speed = np.hypot(q[..., 1], q[..., 2]); peak = np.unravel_index(np.argmax(speed), speed.shape)
        weights = mesh.rf[:, None]*np.ones_like(j)
        weights[-1] *= .5; weights[:, [0, -1]] *= .5
        depth = np.minimum(mesh.radius-mesh.rf[:, None], mesh.half_length-abs(mesh.zf[None, :]))
        samples.append({"time_s": item["time_s"], "maximum_speed_m_s": float(speed.max()),
                        "peak_velocity_location_m": [float(mesh.r[peak[0]]), float(mesh.z[peak[1]])],
                        "kinetic_energy_j": float(np.sum(.5*q[..., 0]*speed**2*mesh.volumes)),
                        "density_min_max_kg_m3": [float(q[..., 0].min()), float(q[..., 0].max())],
                        "current_weighted_rms_depth_m": float(np.sqrt(np.sum(j*j*weights*depth**2)/np.sum(j*j*weights))),
                        "maximum_numerical_eta_over_physical": float(fields[prefix+"_numerical_eta_over_physical"].max())})
    return {"name": row["name"], "method": row["method"], "radial_cells": row["radial_cells"],
            "diagnostics": row["diagnostics"], "samples": samples}


def compare_ends(short, long, window_half_length=.5):
    """Compare identical cells in a central window when actual gas ends move."""
    ar, af = short; br, bf = long
    am, bm = MHDMesh(**ar["fluid_mesh"]), MHDMesh(**br["fluid_mesh"])
    offset = int(round((bm.half_length-am.half_length)/am.dz))
    if am.nr != bm.nr or am.radius != bm.radius or not np.isclose(am.dz, bm.dz) or offset < 0:
        raise ValueError("End comparison requires aligned equal-spacing common cells")
    np.testing.assert_allclose(am.z, bm.z[offset:offset+am.nz], atol=1e-14)
    select = abs(am.z) < window_half_length
    weight = am.volumes[:, select, None]
    rows = []
    for item in ar["fields"]["times"]:
        key = item["prefix"]+"_primitive"
        a = af[key][:, select]; b = bf[key][:, offset:offset+am.nz][:, select]
        bias = ar["source_configuration"]["bias_t"]
        norm = lambda x, y: float(np.sqrt(np.sum((x-y)**2*weight)/np.sum(y*y*weight)))
        rows.append({"time_s": item["time_s"],
                     "velocity_relative_L2": norm(a[..., 1:3], b[..., 1:3]),
                     "magnetic_perturbation_relative_L2": norm(a[..., 4:6]-[0, bias], b[..., 4:6]-[0, bias])})
    return {"short_case": ar["name"], "long_case": br["name"], "window_half_length_m": window_half_length,
            "source_work_change_relative_to_long": ar["diagnostics"]["source_work_j"]/br["diagnostics"]["source_work_j"]-1,
            "history": rows}
