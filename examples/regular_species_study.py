"""Independent rate, conservation and retained transfer-control measurements."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_species import RegularSpecies
from plasma_bench.kernel_induction import KernelInduction
from .regular_species_reference import initialize, exact
from .current_transport_sampling import TransportSampler


def relative(value, target, weight):
    return float(np.sqrt(np.sum(weight*(value-target)**2)/np.sum(weight*target**2)))


def norms(value, target, weight, radius, fluid):
    square = (value-target)**2
    size = target**2
    if value.ndim == weight.ndim+1:
        square, size = square.sum(axis=-1), size.sum(axis=-1)
    rms = np.sqrt(np.sum(weight*size)/np.sum(weight))
    result = dict(relative_L2=float(np.sqrt(np.sum(weight*square)/np.sum(weight*size))),
                  maximum_over_rms=float(np.sqrt(square.max())/rms))
    for label, mask in (('axis_element', radius < fluid.radius/fluid.elements),
                        ('wall_element', radius > fluid.radius*(1-1/fluid.elements))):
        local = weight*mask
        result[label] = float(np.sqrt(np.sum(local*square)/np.sum(local*size)))
    return result


def audit(spec, arrays, order):
    ref, d, _, m, _ = initialize(spec)
    sampler = TransportSampler(d, m.mesh, order)
    s = sampler
    target = exact(ref, s.R, s.Z, d.thermal.gammas, d.thermal.n0, (-m.beta, m.alpha))
    u = arrays['state']
    n = s.scalar(u[..., 0])
    j, jgrad = s.electric(arrays['current'], True)
    results = []
    for k, gamma in enumerate(d.thermal.gammas):
        rate = s.scalar(arrays['fluid_rate'][..., k+4])
        field = s.electric(arrays['pressure_electric'][k])
        entropy = s.scalar(u[..., k+4])/n
        h = E*np.exp(entropy)*(n/d.thermal.n0)**(gamma-1)/(gamma-1)
        pressure_work = float(np.sum(s.weight*np.sum(j*field, axis=-1)))
        thermal_work = float(np.sum(s.weight*h*rate))
        scale = float(np.sum(s.weight*(abs(h*rate)+abs(np.sum(j*field, axis=-1)))))
        results.append(dict(entropy=norms(rate, target['entropy_rate'][k], s.weight, s.R, d),
                            pressure=norms(field, target['pressure'][k], s.weight, s.R, d),
                            continuous_entropy_rate=float(np.sum(s.weight*rate)),
                            continuous_entropy_scale=float(np.sum(s.weight*abs(rate))),
                            pressure_work=pressure_work, continuous_thermal_work=thermal_work,
                            continuous_power_defect=thermal_work-pressure_work,
                            continuous_power_backward=abs(thermal_work-pressure_work)/scale))
    divj = np.trace(jgrad, axis1=-2, axis2=-1)
    charge = float(np.sqrt(ref.radius**2*np.sum(s.weight*divj**2)
                          /np.sum(s.weight*np.sum(j*j, axis=-1))))
    return dict(species=results, strong_charge_divergence_scaled=charge)


def run_case(spec):
    ref, d, u, m, field = initialize(spec)
    g = m.geometry
    current = field.current(m)
    species = RegularSpecies(m)
    pair = species.pair(u, current)
    target = exact(ref, d.r, d.z, d.thermal.gammas, d.thermal.n0, (-m.beta, m.alpha))
    electric = np.stack([np.concatenate([solve(v) for solve, v in
                         zip(g.solves, g.split_electric(force))])
                         for force in pair['species_forcing']])
    f = d.fields(u)
    records = []
    for k in range(2):
        rate = pair['fluid_rate'][..., k+4]
        work = float(np.sum(d.volume*f['specific_internal'][..., k]*rate))
        current_work = float(current@pair['species_forcing'][k])
        scale = float(np.sum(d.volume*abs(f['specific_internal'][..., k]*rate))
                      + abs(current)@abs(pair['species_forcing'][k]))
        entropy_scale = float(np.sum(d.volume*abs(rate)))
        quadratic_scale = float(np.sum(d.volume*abs(f['entropy'][..., k]*rate)))
        records.append(dict(point_rate=norms(rate, target['entropy_rate'][k], d.volume, d.r, d),
                            entropy_backward=float(abs(np.sum(d.volume*rate))/entropy_scale),
                            quadratic_backward=float(abs(np.sum(d.volume*f['entropy'][..., k]*rate))/quadratic_scale),
                            thermal_work=work, current_work=current_work,
                            work_backward=abs(work-current_work)/scale))
    solved = KernelInduction(m, tolerance=1e-10, max_iterations=3000).solve(pair['reduced_forcing'])
    full_electric = solved.electric+pair['potential_electric']
    original_ohm = m.mass@full_electric-m.inertia@solved.current-pair['forcing']
    original_scale = abs(m.mass@full_electric)+abs(m.inertia@solved.current)+abs(pair['forcing'])
    nt = np.zeros_like(u[..., 0])
    thermal = sum(r['thermal_work'] for r in records)
    magnetic = float(np.sum(m.mesh.volumes*m.magnetic_energy_rate(field, solved.magnetic)))
    kinetic = float(np.sum(m.mesh.volumes*m.current_energy_rate(current, solved.current, nt)['total']))
    point = RegularSpecies(m, transfer='point')
    controls = {}
    for label, operator, skew in (('point_unskew', point, False), ('point_skew', point, True)):
        rate = m.beta/E*operator.transport(current, f['entropy'][..., 0], skew=skew)/d.volume
        controls[label] = relative(rate, target['entropy_rate'][0], d.volume)
    arrays = dict(state=u, current=current, magnetic=field.packed_increment(g),
                  fluid_rate=pair['fluid_rate'], pressure_electric=electric,
                  species_forcing=pair['species_forcing'], magnetic_rate=solved.magnetic,
                  current_rate=solved.current, electric=full_electric,
                  reduced_electric=solved.electric, potential_electric=pair['potential_electric'])
    record = dict(spec=spec, status='complete', species=records, controls=controls,
                  thermal_work=thermal, magnetic_work=magnetic, current_kinetic_work=kinetic,
                  total_grid_power=thermal+magnetic+kinetic,
                  total_grid_power_backward=abs(thermal+magnetic+kinetic)/(abs(thermal)+abs(magnetic)+abs(kinetic)),
                  induction_residual=solved.residual, residual_components=solved.residual_components)
    record['original_ohm_backward'] = float(np.linalg.norm(original_ohm)/np.linalg.norm(original_scale))
    return record, arrays
