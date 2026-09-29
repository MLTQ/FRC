"""Physical nested-mesh and same-grid temporal differences for warm trajectories."""
import numpy as np
from .warm_hall_norms import sample_q1, pressure_difference
from .compatible_hall_audit import magnetic_product


def species_pressure(record, arrays, prefix):
    if record['configuration']['kind'] != 'pressure_evolve':
        return arrays[prefix+'_electron_pressure']
    gamma = record['configuration']['reference']['gamma']
    return (arrays['number_density']/1e20)**(gamma-1)*arrays[prefix+'_entropy']


def magnetic_fields(arrays, prefix, mesh):
    if prefix+'_conserved' in arrays:
        return arrays[prefix+'_psi'], arrays[prefix+'_conserved'][..., 5]
    return np.zeros((mesh.nr+1, mesh.nz+1)), arrays[prefix+'_g']


def difference(records, arrays, meshes, coarse_name, fine_name):
    c, f = meshes[coarse_name], meshes[fine_name]
    assert c.radius == f.radius and c.half_length == f.half_length
    assert f.nr % c.nr == 0 and f.nz % c.nz == 0
    np.testing.assert_allclose(c.rf, f.rf[::f.nr//c.nr], rtol=2e-15, atol=1e-18)
    assert abs(records[coarse_name]['end_time_s']/records[fine_name]['end_time_s']-1) < 2e-13
    a, b = arrays[coarse_name], arrays[fine_name]
    pe = pressure_difference(c, species_pressure(records[coarse_name], a, 'final'), f,
                             species_pressure(records[fine_name], b, 'final'), species_pressure(records[fine_name], b, 'initial'))
    cp, cg = magnetic_fields(a, 'final', c)
    fp, fg = magnetic_fields(b, 'final', f)
    ip, ig = magnetic_fields(b, 'initial', f)
    # Q1 flux and P0 g prolong exactly into the nested fine representation.
    prolonged = sample_q1(c, cp[:, :-1], .5*f.rf[:, None]**2, f.zf[None])
    prolonged[:, -1] = prolonged[:, 0]
    cg = np.repeat(np.repeat(cg, f.nr//c.nr, axis=0), f.nz//c.nz, axis=1)
    dp, dg = prolonged-fp, cg-fg
    error = np.sum(f.volumes*magnetic_product(f, dp, dg, dp, dg))
    change = np.sum(f.volumes*magnetic_product(f, fp-ip, fg-ig, fp-ip, fg-ig))
    return dict(coarse=coarse_name, reference=fine_name, electron_pressure=pe,
                magnetic_relative_to_fine_field_change=float(np.sqrt(error/change)) if change else None,
                magnetic_energy_norm_difference_j=float(.5*error))


def comparisons(records, arrays, meshes):
    temporal = {}
    for label, counts, finest in (('pressure', (10, 20, 40), 320), ('wave', (8, 16, 32), 256), ('nonlinear', (10, 20, 40), 160)):
        values = [difference(records, arrays, meshes, f'{label}_time_{n}', f'{label}_time_{finest}') for n in counts]
        p = [v['electron_pressure']['rms_difference_pa'] for v in values]
        b = [v['magnetic_relative_to_fine_field_change'] for v in values]
        temporal[label] = dict(comparisons=values, pressure_difference_ratios=[p[0]/p[1], p[1]/p[2]],
                               magnetic_difference_ratios=[b[0]/b[1], b[1]/b[2]])
    spatial = {label: [difference(records, arrays, meshes, f'{label}_{n}', f'{label}_{fine}') for n in counts]
               for label, counts, fine in (('pressure_evolve', (8, 16, 32), 64), ('nonlinear', (8, 16), 32))}
    controls = dict(amplitudes={name: r['metrics'] for name, r in records.items() if name.startswith('wave_amplitude_')},
                    reverse_waves={f'branch_{k}': {key: records[f'wave_reverse_b{k}']['metrics'][key]-records[f'wave_b{k}_16']['metrics'][key]
                                                  for key in records[f'wave_b{k}_16']['metrics']} for k in range(3)})
    a, b = arrays['wave_b2_16']['final_conserved'], arrays['wave_cold_limit']['final_conserved']
    controls['near_cold_electrons'] = dict(max_absolute_conserved_difference_by_component=np.max(abs(a-b), axis=(0, 1)).tolist(),
                                           electron_fraction=records['wave_cold_limit']['configuration']['reference']['electron_fraction'])
    a, b = arrays['trace_diffusion_128'], arrays['trace_diffusion_half_time']
    weight = meshes['trace_diffusion_128'].volumes
    controls['diffusion_two_half_steps'] = dict(g_relative_change=float(np.sqrt(np.sum(weight*(a['final_g']-b['final_g'])**2)/np.sum(weight*b['final_g']**2))),
                                                one_step_metrics=records['trace_diffusion_128']['metrics'], two_step_metrics=records['trace_diffusion_half_time']['metrics'])
    return dict(time=temporal, mesh=spatial, controls=controls)
