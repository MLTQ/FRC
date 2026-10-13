"""Reconstruct actual replayed trajectory states for an honest 3D display."""
import argparse
import base64
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from plasma_bench.regular_current import RegularCurrent
from .current_transport_sampling import bases
from .split_mode_audit import fixture
from .trajectory_tangent_probe import verified_inputs
from .implicit_plasma_records import provenance, freeze, sha


class FieldSampler:
    def __init__(self, mesh, r, z):
        self.mesh = mesh
        self.r = np.asarray(r)
        self.nx, self.dx, _, self.ex, _ = bases(.5*mesh.rf**2, .5*self.r**2, 5)
        self.nz, self.dz, _, self.ez, _ = bases(mesh.zf, np.asarray(z), 5, True)

    def fields(self, magnetic, current, bias):
        m = self.mesh
        nodes = (m.nr+1)*m.nz
        p = magnetic[:nodes].reshape(m.nr+1, m.nz)
        g = magnetic[nodes:].reshape(m.nr, m.nz)
        r = self.r[:, None]
        b = np.stack((-(self.nx@p@self.dz.T)/r,
                      r*(self.ex@g@self.ez.T), self.dx@p@self.nz.T+bias), axis=-1)
        count, nt = m.nr*m.nz, (m.nr-1)*m.nz
        jr = current[:count].reshape(m.nr, m.nz)
        jt = np.pad(current[count:count+nt].reshape(m.nr-1, m.nz), ((1, 1), (0, 0)))
        jz = current[count+nt:].reshape(m.nr+1, m.nz)
        j = np.stack((r*(self.ex@jr@self.nz.T),
                      (self.nx@(m.rf[:, None]*jt)@self.nz.T)/r,
                      self.nx@jz@self.ez.T), axis=-1)
        if not np.isfinite(b).all() or not np.isfinite(j).all():
            raise ValueError('Finite reconstructed vectors required')
        return b, j


def errors(value, target):
    return dict(relative_L2=float(np.linalg.norm(value-target)/np.linalg.norm(target)),
                max_absolute=float(np.max(abs(value-target))))


def trace(r, z, magnetic, radius, scale, steps):
    if np.min(magnetic[..., 2]) <= 0:
        raise ValueError('This z-parameterized display requires strictly positive Bz')
    interp = RegularGridInterpolator((r, z), magnetic, bounds_error=True)
    seeds = np.array([(radius*rad, phi) for rad in (.12, .28, .44, .60, .76, .91)
                      for phi in np.arange(8)*2*np.pi/8])
    axial = np.linspace(z[0], z[-1], steps+1)
    path = np.empty((steps+1, len(seeds), 3))
    strengths = np.empty((steps+1, len(seeds)))
    q = seeds.copy()
    def slope(value, zz):
        b = interp(np.column_stack((value[:, 0], np.full(len(value), zz))))
        if np.any(b[:, 2] <= 0): raise ValueError('Nonpositive axial field on a trace')
        return np.column_stack((b[:, 0]/b[:, 2], b[:, 1]/(value[:, 0]*b[:, 2])))
    for k, zz in enumerate(axial):
        rr, angle = q.T
        path[k] = np.column_stack((np.full(len(q), zz), rr*np.cos(angle), rr*np.sin(angle)))*scale
        strengths[k] = np.linalg.norm(interp(np.column_stack((rr, np.full(len(q), zz)))), axis=-1)*1000
        if k == steps: break
        dz = axial[k+1]-zz
        k1 = slope(q, zz)
        k2 = slope(q+.5*dz*k1, zz+.5*dz)
        k3 = slope(q+.5*dz*k2, zz+.5*dz)
        k4 = slope(q+dz*k3, axial[k+1])
        q += dz*(k1+2*k2+2*k3+k4)/6
    return path, strengths


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Fresh display data destination required')
    trajectory = Path('outputs/trajectory_refined_t60_minus')
    verification = Path('outputs/trajectory_refined_t60_minus_replay_verification')
    original, inputs = verified_inputs(trajectory, verification)
    if original['completed_steps'] != 17 or original['success'] or original['failed_step'] != 18:
        raise ValueError('Exact documented terminal trajectory required')
    identity = provenance(['examples.trajectory_display_data'])
    destination.mkdir(parents=True)
    freeze(destination, identity)
    mode, d, _, system, bias, _, _ = fixture(original['spec'])
    if mode.input_sha256 != original['input_sha256']: raise ValueError('Original fixture required')
    mesh, geo = system.geometry.mesh, system.geometry
    nr, nz = 256, 512
    r = (np.arange(nr)+.5)*d.radius/nr
    z = ((np.arange(nz)+.5)/nz-.5)*d.length
    sample = FieldSampler(mesh, r, z)
    production = FieldSampler(mesh, geo.radial_points, geo.axial_points)
    check_r = r[1:-1:2]+.25*d.radius/nr
    check_z = z[1:-1:2]+.25*d.length/nz
    check = FieldSampler(mesh, check_r, check_z)
    cr, cz = np.meshgrid(check_r, check_z, indexing='ij')
    check_points = np.column_stack((cr.ravel(), cz.ravel()))
    scale = 1/np.hypot(d.radius, .5*d.length)
    frames, measures = [], []
    for step in (0, 10, 15, 17):
        record = original['records'][step-1] if step else None
        name = record['arrays_file'] if step else 'initial.npz'
        with np.load(trajectory/name) as raw:
            state, magnetic = raw['state'].copy(), raw['magnetic'].copy()
        model = RegularCurrent(system.mapping, state[..., 0], system.ion_mass)
        current = model.current(magnetic)
        b, j = sample.fields(magnetic, current, bias.axial_bias)
        qb, qj = production.fields(magnetic, current, bias.axial_bias)
        pb = geo.magnetic_samples(*geo.unpack(magnetic)); pb[..., 2] += bias.axial_bias
        pj = geo.sample_electric(current)
        unpack = lambda v: np.stack([geo._unlayout(v[..., k]) for k in range(3)], axis=-1)
        agreement = dict(B=errors(qb, unpack(pb)), J=errors(qj, unpack(pj)))
        if max(v['relative_L2'] for v in agreement.values()) > 1e-10:
            raise ValueError('Independent field reconstruction disagrees with production')
        cb, cj = check.fields(magnetic, current, bias.axial_bias)
        interpolation = {key: errors(RegularGridInterpolator((r, z), value)(check_points).reshape(target.shape), target)
                         for key, value, target in [('B', b, cb), ('J', j, cj)]}
        points, strength = trace(r, z, b, d.radius, scale, 400)
        fine_points, _ = trace(r, z, b, d.radius, scale, 800)
        trace_error = float(np.max(np.linalg.norm(points-fine_points[::2], axis=-1))/scale)
        jm, bm = np.linalg.norm(j, axis=-1), np.linalg.norm(b, axis=-1)
        texture = np.stack((jm, bm, np.zeros_like(jm), b[..., 0]), axis=-1).transpose(1, 0, 2).astype('<f4')
        payload = dict(step=step, time_us=(record['time'] if step else 0.)*1e6,
                       bias_mT=bias.axial_bias*1000, sampled_Jmax_A_m2=float(jm.max()),
                       sampled_Bmin_mT=float(bm.min()*1000), sampled_Bmax_mT=float(bm.max()*1000),
                       sampled_Bzmin_mT=float(b[..., 2].min()*1000))
        frames.append(dict(**payload, texture=base64.b64encode(texture.tobytes()).decode(),
                           lines3d=[dict(points=points[:, k].tolist(), b=strength[:, k].tolist()) for k in range(points.shape[1])]))
        file = destination/f'step_{step:04}.npz'
        np.savez_compressed(file, r=r, z=z, B=b, J=j, current_coefficients=current,
                            line_points=points, line_strength_mT=strength)
        measures.append(dict(**payload, source_file=str(trajectory/name), source_sha256=sha(trajectory/name),
                             arrays_file=file.name, arrays_sha256=sha(file),
                             independent_vs_production=agreement, display_interpolation=interpolation,
                             rk4_half_step_max_position_difference_m=trace_error))
        print('saved field display', step, payload['time_us'], payload['sampled_Jmax_A_m2'], flush=True)
    data = dict(frames=frames, texture_size=[nr, nz], extent=[d.radius*scale, .5*d.length*scale],
                uniform_display_scale_per_m=scale, physical_radius_m=d.radius, physical_length_m=d.length,
                jmax=max(v['sampled_Jmax_A_m2'] for v in measures),
                bmin_mT=min(v['sampled_Bmin_mT'] for v in measures), bmax_mT=max(v['sampled_Bmax_mT'] for v in measures),
                trajectory=str(trajectory), trajectory_source_sha256=original['source_sha256'],
                interpretation='Axisymmetric contact verification fixture, revolved for display. Step 18 rejected; no FRC formation or 3D stability claim.')
    path = destination/'data.json'; path.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    if any(sha(p) != h for p, h in inputs.items()) or provenance(['examples.trajectory_display_data']) != identity:
        raise ValueError('Original inputs or display source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], input_sha256=inputs,
                  data_sha256=sha(path), frames=measures, spec=original['spec'],
                  scope=data['interpretation'], source_trajectory_terminal=True, full_duration_success=False,
                  failed_step=18, requested_duration_us=100., sampled_time_max_us=original['final_time']*1e6,
                  sampling='Independent degree-five magnetic/current basis; linear display interpolation; no time interpolation.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    run(parser.parse_args().destination)
