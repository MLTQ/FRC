"""Resolution-limited critical points, closed contours and separatrix estimates."""
import numpy as np
import contourpy
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import root
from matplotlib.path import Path
from .numerics import field_from_flux


def _volume(path):
    x, y = path[:, 0], path[:, 1]  # z, r
    cross = x*np.roll(y, -1)-np.roll(x, -1)*y
    return float(abs(2*np.pi*np.sum((y+np.roll(y, -1))*cross)/6))


def analyze_flux(grid, psi, region=None):
    """Topology inside region (usually plasma), independent of field reversal."""
    r, z, p = grid.r, grid.z, np.asarray(psi)
    if region is not None:
        ri = np.flatnonzero(region.contains(r, 0))
        zi = np.flatnonzero(region.contains(0, z))
        r, z, p = r[ri], z[zi], p[np.ix_(ri, zi)]
    result = {"o_points": [], "x_points": [], "axis_nulls": [],
              "closed_flux": False, "closed_contour_count": 0,
              "largest_closed_radius_m": None, "largest_closed_volume_m3": None,
              "separatrix_radius_m": None, "separatrix_volume_m3": None,
              "separatrix_method": None,
              "resolution_note": "Sampled contours and interpolated critical points; unresolved features may be missed."}
    scale = float(np.ptp(p))
    if min(p.shape) < 5 or scale < 1e-25:
        return result
    f = RectBivariateSpline(r, z, p/scale)
    gr, gz = np.gradient(p/scale, r, z, edge_order=2)
    mins_r = np.minimum.reduce([gr[:-1, :-1], gr[1:, :-1], gr[:-1, 1:], gr[1:, 1:]])
    maxs_r = np.maximum.reduce([gr[:-1, :-1], gr[1:, :-1], gr[:-1, 1:], gr[1:, 1:]])
    mins_z = np.minimum.reduce([gz[:-1, :-1], gz[1:, :-1], gz[:-1, 1:], gz[1:, 1:]])
    maxs_z = np.maximum.reduce([gz[:-1, :-1], gz[1:, :-1], gz[:-1, 1:], gz[1:, 1:]])
    candidates = np.argwhere((mins_r <= 0) & (maxs_r >= 0) & (mins_z <= 0) & (maxs_z >= 0))
    locations = []
    for i, j in candidates:
        if i == 0:
            continue
        solution = root(lambda q: [f.ev(*q, dx=1), f.ev(*q, dy=1)], [(r[i]+r[i+1])/2, (z[j]+z[j+1])/2])
        rr, zz = solution.x
        if not solution.success or not (r[1] < rr < r[-2] and z[1] < zz < z[-2]):
            continue
        if np.linalg.norm(solution.fun)*max(r[-1], z[-1]-z[0]) > 1e-6:
            continue
        if any(np.hypot((rr-a)/grid.dr, (zz-b)/grid.dz) < .5 for a, b in locations):
            continue
        locations.append((rr, zz))
        hessian = np.array([[f.ev(rr, zz, dx=2), f.ev(rr, zz, dx=1, dy=1)],
                            [f.ev(rr, zz, dx=1, dy=1), f.ev(rr, zz, dy=2)]])
        eig = np.linalg.eigvalsh(hessian)
        if np.min(np.abs(eig)) < 1e-9*max(np.max(np.abs(eig)), 1):
            continue  # Degenerate critical points are not confidently classified.
        key = "o_points" if eig[0]*eig[1] > 0 else "x_points"
        result[key].append({"r_m": float(rr), "z_m": float(zz), "psi_wb": float(f.ev(rr, zz)*scale)})
    # Flux gradient vanishes identically on axis; classify physical Bz zeros instead.
    _, full_bz = field_from_flux(grid, psi)
    axis = full_bz[0] if region is None else full_bz[0, zi]
    for j in range(len(z)-1):
        if axis[j]*axis[j+1] < 0 or (axis[j] == 0 and axis[j+1] != 0):
            zz = z[j]-axis[j]*(z[j+1]-z[j])/(axis[j+1]-axis[j])
            result["axis_nulls"].append({"r_m": 0.0, "z_m": float(zz), "psi_wb": 0.0})
    cg = contourpy.contour_generator(x=z, y=r, z=p, line_type="Separate")
    closed = []
    levels = list(np.linspace(float(p.min()), float(p.max()), 34)[1:-1])
    for o in result["o_points"]:
        levels.extend(o["psi_wb"]*fraction for fraction in (.1, .5, .9))
    for level in levels:
        for line in cg.lines(level):
            if len(line) < 8 or not np.allclose(line[0], line[-1], rtol=0, atol=1e-10):
                continue
            if not any(Path(line).contains_point((o["z_m"], o["r_m"])) for o in result["o_points"]):
                continue
            closed.append(line)
    result["closed_contour_count"] = len(closed)
    result["closed_flux"] = bool(closed)
    if closed:
        largest = max(closed, key=_volume)
        result["largest_closed_radius_m"] = float(largest[:, 1].max())
        result["largest_closed_volume_m3"] = _volume(largest)
    # Axis-attached separatrix: psi=0 arc connecting two resolved axis nulls.
    # Remove the identically-zero axis row before contour extraction.
    if closed and len(result["axis_nulls"]) >= 2:
        cg_axis = contourpy.contour_generator(x=z, y=r[1:], z=p[1:], line_type="Separate")
        separatrices = []
        for arc in cg_axis.lines(0):
            if len(arc) < 6 or max(arc[0, 1], arc[-1, 1]) > r[1]+1e-8:
                continue
            zends = sorted([arc[0, 0], arc[-1, 0]])
            nulls = [a["z_m"] for a in result["axis_nulls"]]
            if any(min(abs(end-n) for n in nulls) > 3*grid.dz for end in zends):
                continue
            polygon = np.vstack((arc, [arc[-1, 0], 0], [arc[0, 0], 0], arc[0]))
            if any(Path(polygon).contains_point((o["z_m"], o["r_m"])) for o in result["o_points"]):
                separatrices.append(polygon)
        if separatrices:
            sep = max(separatrices, key=_volume)
            result.update(separatrix_radius_m=float(sep[:, 1].max()),
                          separatrix_volume_m3=_volume(sep),
                          separatrix_method="axis-flux contour closed to axis; grid estimate")
    return result
