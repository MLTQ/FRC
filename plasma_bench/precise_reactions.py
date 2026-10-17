"""Conservative reaction gather/return with distinct precise material supports."""
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from .mass_scaling import solve_mass_banded
from .projection_energy import projected_frame_energy
from .bounded_reaction_projection import bounded_moments
from .constrained_velocity_projection import constrained_velocity
from .variance_projection import variance_projection
from .velocity_extension import inherited_velocity_faces
from .precise_geometry import decimal, widths, fraction, quadrature, interpolate, geometry_precision
from .precise_material import PreciseMaterialSubcells, PreciseMaterialState, PreciseMaterialSet
from .material_species_set import MaterialSetState
from .material_reactions import MaterialReactionOverlap, ReactionOverlapState
from .reactive_material import ReactiveMaterial, MaterialReservoirs


class PreciseReactionOverlap(MaterialReactionOverlap):
    @geometry_precision
    def from_material(self, model, state, species):
        species = np.asarray(species)
        if (species.shape != (model.species,) or not np.issubdtype(species.dtype, np.integer)
                or np.any(species < 0) or np.any(species >= len(self.number))):
            raise ValueError('Atomic species index for every precise support required')
        if not np.array_equal(model.mass, self.source.mass[species]) or not np.array_equal(model.charge, self.charge[species]):
            raise ValueError('Identical material/reaction species masses and charges required')
        model.validate(state)
        edges = np.unique(np.r_[model.bounds, state.faces])
        points, measure = quadrature(edges, self.quadrature)
        sample = model.sample(state, points.ravel())
        rows, network = [], self.source.network
        for k in range(points.size):
            atomic = np.zeros(self.source.base_size)
            temperature = np.zeros(network.count+1)
            velocity = np.zeros((network.count+1, 3))
            active = sample['number'][:, k] > 0
            ids = species[active]
            if len(set(ids)) != len(ids):
                raise ValueError('Overlapping same-species blocks need an explicit mixing closure')
            atomic[self.number[ids]] = sample['number'][active, k]/self.density
            if np.any(atomic[self.number[ids]] == 0):
                raise ValueError('Positive density below atomic normalization range')
            temperature[ids] = sample['temperature_ev'][active, k]
            velocity[ids] = sample['velocity'][active, k]
            atomic[network.thermal] = 1.5*atomic[network.electrons]*temperature[-1]
            rows.append(self.source.initial(atomic, temperature[:-1], velocity))
        return ReactionOverlapState(edges, points, measure,
            np.asarray(rows).reshape(points.shape+(self.source.size,)), model.left_electric)

    @geometry_precision
    def totals(self, state):
        f, source, network = self.fields(state), self.source, self.source.network
        sums = self.density*np.sum(state.weights[..., None]*state.moments, axis=(0, 1))
        numbers = sums[self.number]
        relative_p = source.momentum_scale*sums[source.momentum].sum(axis=0)
        mass = source.mass@numbers
        energy = EV*(sums[source.energy].sum()+network.energy_ev@sums[:network.count]+
                     sums[network.photons]-sums[network.heating])
        energy += self.frame@relative_p+.5*mass*(self.frame@self.frame)
        charge = f['charge_c_m3'].reshape(state.points.shape)
        electric = np.r_[state.left_electric, state.left_electric+
                         np.cumsum(np.sum(state.weights*charge, axis=1))/EPS0]
        field = EPS0/6*np.sum(widths(state.edges)*(electric[:-1]**2+
                             electric[:-1]*electric[1:]+electric[1:]**2))
        return dict(numbers_m2=numbers, nuclei_m2=float(numbers[:-1].sum()),
            charge_c_m2=float(EV*(numbers@self.charge)), mass_kg_m2=float(mass),
            momentum_kg_per_m_s=relative_p+mass*self.frame,
            matter_binding_photon_energy_j_m2=float(energy), electric_energy_j_m2=float(field),
            energy_j_m2=float(energy+field), photon_energy_j_m2=float(EV*sums[network.photons]),
            electric_edges_v_m=electric)


@geometry_precision
def project_precise_reactions(overlap, state, *, velocity_cells, coalesce_ulps=0,
        velocity_knots=None, energy_projection='local_total', projection_method='galerkin', extend_velocity=False, cfl=.15, maximum_plasma_phase=.04,
        numerical_viscosity=(.5, 2.), viscosity_limiter='minmod', newborn_spacing_fraction=0.):
    if not isinstance(coalesce_ulps, int) or coalesce_ulps < 0:
        raise ValueError('Nonnegative integer explicitly declared ULP cleanup required')
    if energy_projection not in ('local_total', 'projected_frame'):
        raise ValueError('Explicit local_total or projected_frame energy projection required')
    if projection_method not in ('galerkin', 'bounded_moments', 'constrained_velocity', 'variance'):
        raise ValueError('Explicit supported reaction projection method required')
    if projection_method != 'galerkin' and energy_projection != 'local_total':
        raise ValueError('Alternative momentum methods require their paired energy transfer')
    if extend_velocity and projection_method == 'bounded_moments':
        raise ValueError('Convex bounded transfer cannot use extrapolating velocity shapes')
    source, network, values = overlap.source, overlap.source.network, state.moments
    present = np.any(values[..., overlap.number] > 0, axis=1)
    protected = np.r_[True, np.any(present[:-1] != present[1:], axis=1), True]
    if velocity_knots:
        protected |= np.isin(state.edges, np.concatenate(list(velocity_knots.values())))
    # Keep the earlier declared float-coordinate cleanup scale, not a new
    # physical-width averaging option. Every species endpoint stays protected.
    tolerance = decimal(coalesce_ulps*np.spacing(max(abs(float(x)) for x in state.edges)))
    kept, i = list(range(len(state.edges))), 0
    while tolerance and i < len(kept)-1:
        a, b = kept[i:i+2]
        if state.edges[b]-state.edges[a] <= tolerance:
            if not protected[b]:
                kept.pop(i+1)
                continue
            if not protected[a]:
                kept.pop(i)
                i = max(0, i-1)
                continue
        i += 1
    edges = state.edges[kept]
    groups = np.searchsorted(kept, np.arange(len(state.edges)-1), side='right')-1
    right = fraction(state.points, edges[groups, None], edges[groups+1, None])
    count, models, states, ids, records = len(edges)-1, [], [], [], []
    thermal_gain, energy_l1, maximum_temperature_change = 0., 0., 0.
    for species, mass in enumerate(source.mass):
        scale = float(values[..., overlap.number[species]].max())
        if scale == 0:
            continue
        population_scale = overlap.density*scale
        nn, ee, source_thermal = np.zeros(count), np.zeros(count), np.zeros(count)
        ll, lr = np.zeros((count, 3)), np.zeros((count, 3))
        el, er = np.zeros(count), np.zeros(count)
        np.add.at(nn, groups, np.sum(state.weights*(values[..., overlap.number[species]]/scale), axis=1))
        pm = values[..., source.momentum[species]]/scale*source.momentum_scale
        em = values[..., source.energy[species]]/scale*EV
        np.add.at(ll, groups, np.sum(state.weights[..., None]*pm*(1-right[..., None]), axis=1))
        np.add.at(lr, groups, np.sum(state.weights[..., None]*pm*right[..., None], axis=1))
        np.add.at(ee, groups, np.sum(state.weights*em, axis=1))
        if projection_method == 'bounded_moments':
            np.add.at(el, groups, np.sum(state.weights*em*(1-right), axis=1))
            np.add.at(er, groups, np.sum(state.weights*em*right, axis=1))
        number = values[..., overlap.number[species]]
        source_velocity = np.divide(values[..., source.momentum[species]], number[..., None],
            out=np.zeros_like(pm), where=number[..., None] > 0)*(source.momentum_scale/mass)
        source_kinetic = (number/scale)*(.5*mass*np.sum(source_velocity*source_velocity, axis=-1))
        np.add.at(source_thermal, groups, np.sum(state.weights*(em-source_kinetic), axis=1))
        occupied = nn > 0
        expected = np.zeros(count, bool)
        np.logical_or.at(expected, groups, present[:, species])
        if np.any(expected & ~occupied):
            raise ValueError('Positive normalized amount is unrepresentable; no discarded species')
        pairs = np.flatnonzero(np.diff(np.r_[False, occupied, False])).reshape(-1, 2)
        for first, last in pairs:
            faces = edges[first:last+1]
            inherited = None
            extended = False
            if velocity_knots is not None and species in velocity_knots:
                old = velocity_knots[species]
                inherited, extended = inherited_velocity_faces(faces, old, extend=extend_velocity)
            model = PreciseMaterialSubcells(nn[first:last], mass, overlap.charge[species],
                (state.edges[0], state.edges[-1]), faces, velocity_cells, number_scale=population_scale,
                velocity_faces=inherited, extend_velocity=extended, newborn_spacing_fraction=newborn_spacing_fraction,
                left_electric=state.left_electric, cfl=cfl, maximum_plasma_phase=maximum_plasma_phase,
                numerical_viscosity=numerical_viscosity, viscosity_limiter=viscosity_limiter)
            load = np.zeros((len(faces), 3))
            load[:-1] += ll[first:last]
            load[1:] += lr[first:last]
            transfer_record = {}
            if projection_method == 'bounded_moments':
                velocity, thermal, redistribution, transfer_record = bounded_moments(
                    model, ll[first:last], lr[first:last], el[first:last], er[first:last])
            elif projection_method == 'variance':
                velocity, thermal, redistribution, transfer_record = variance_projection(model,
                    state.points.ravel(), state.weights.ravel(), source_velocity.reshape(-1, 3),
                    (number/scale).ravel(), (em-source_kinetic).ravel(),
                    np.broadcast_to(groups[:, None]-first, state.points.shape).ravel(),
                    load.sum(axis=0), ee[first:last])
            elif projection_method == 'constrained_velocity':
                velocity, transfer_record = constrained_velocity(model,
                    state.points.ravel(), state.weights.ravel(), source_velocity.reshape(-1, 3),
                    (number > 0).ravel(), load.sum(axis=0), ee[first:last])
            else:
                velocity = solve_mass_banded(model.matrix[0], model.restrict_load(load))
            vl, vr = model.fine_velocity(velocity[None, ...])
            kinetic = mass*model.relative_numbers/6*np.sum(vl*vl+vl*vr+vr*vr, axis=1)
            if projection_method in ('bounded_moments', 'variance'):
                pass
            elif energy_projection == 'projected_frame':
                thermal, redistribution = projected_frame_energy(ee[first:last], ll[first:last],
                    lr[first:last], vl, vr, kinetic)
            else:
                thermal, redistribution = ee[first:last]-kinetic, np.zeros(last-first)
            result = PreciseMaterialState(faces[None, :].copy(), velocity[None, ...],
                                          thermal[None, :])
            model.validate(result)
            gain = thermal-source_thermal[first:last]
            transfer_record.update(
                maximum_temperature_change_ev=float(abs((model.gamma-1)*(gain/model.relative_numbers)/EV).max()),
                minimum_temperature_ev=float(((model.gamma-1)*thermal/model.relative_numbers/EV).min()))
            thermal_gain += float(gain.sum()*population_scale)
            energy_l1 += float(abs(redistribution).sum()*population_scale)
            maximum_temperature_change = max(maximum_temperature_change,
                float(abs((model.gamma-1)*(gain/model.relative_numbers)/EV).max()))
            records.append(dict(species=species, cells=last-first, minimum_width_m=float(widths(faces).min()),
                velocity_cells=model.velocity_cells, inherited_velocity_faces=inherited is not None,
                extended_velocity_support=extended,
                support_decimal_m=[str(faces[0]), str(faces[-1])], **transfer_record))
            models.append(model)
            states.append(result)
            ids.append(species)
    if not models:
        raise ValueError('No occupied support; retain exact vacuum without a matter solver')
    model = PreciseMaterialSet(models, ids, reference_velocity=overlap.frame)
    result = MaterialSetState(tuple(states))
    model.validate(result)
    amounts = overlap.density*np.sum(state.weights[..., None]*values[..., overlap.number], axis=1)
    old = np.r_[state.left_electric, state.left_electric+np.cumsum(EV*(amounts@overlap.charge))/EPS0]
    new_edges, new, _, new_energy = model.field(result)
    delta = interpolate(state.edges, new_edges, new)-old
    h = widths(state.edges)
    old_energy = EPS0/6*np.sum(h*(old[:-1]**2+old[:-1]*old[1:]+old[1:]**2))
    record = dict(blocks=records, removed_roundoff_edges=len(state.edges)-len(edges),
        projection_method=projection_method,
        energy_projection=energy_projection, projection_thermal_gain_j_m2=thermal_gain,
        projection_energy_redistribution_l1_j_m2=energy_l1,
        maximum_projection_temperature_change_ev=maximum_temperature_change,
        geometric_tolerance_m=float(tolerance), removed_density_edges=0, minimum_density_width_m=0.,
        projection_field_energy_change_j_m2=float(new_energy-old_energy),
        projection_field_maximum_change_v_m=float(abs(delta).max()),
        projection_field_l2_change_v_sqrt_m=float(np.sqrt(np.sum(h*(delta[:-1]**2+delta[:-1]*delta[1:]+delta[1:]**2))/3)),
        photon_energy_j_m2=float(EV*overlap.density*np.sum(state.weights*values[..., network.photons])),
        input_heat_j_m2=float(EV*overlap.density*np.sum(state.weights*values[..., network.heating])),
        binding_energy_j_m2=float(EV*np.sum(amounts[:, :network.count]*network.energy_ev)))
    return model, result, record


class PreciseReactiveMaterial(ReactiveMaterial):
    def __init__(self, source, reference_density, *, quadrature=3, preserve_velocity_knots=True, **controls):
        super().__init__(source, reference_density, quadrature=quadrature, **controls)
        self.overlap = PreciseReactionOverlap(source, reference_density, quadrature=quadrature)
        self.preserve_velocity_knots = bool(preserve_velocity_knots)
        self.last_velocity_knots = None

    def source_step(self, model, state, reservoirs, dt):
        self.last_overlap = None
        self.last_velocity_knots = None
        gathered = self.overlap.from_material(model, state, model.species_indices)
        self.last_overlap = self.overlap.step(gathered, dt)
        knots = None
        if self.preserve_velocity_knots:
            knots = {}
            for index, block, values in zip(model.species_indices, model.models, state.blocks):
                knots.setdefault(int(index), []).append(values.faces[0, block.velocity_indices])
            knots = {index: np.unique(np.concatenate(values)) for index, values in knots.items()}
        self.last_velocity_knots = knots
        result, advanced, record = project_precise_reactions(self.overlap, self.last_overlap,
            velocity_knots=knots, **self.controls)
        self.last_projection = record
        ledger = MaterialReservoirs(reservoirs.photons_j_m2+record['photon_energy_j_m2'],
            reservoirs.input_heat_j_m2+record['input_heat_j_m2'],
            reservoirs.projection_field_change_j_m2+record['projection_field_energy_change_j_m2'],
            reservoirs.projection_thermal_gain_j_m2+record['projection_thermal_gain_j_m2'],
            reservoirs.projection_energy_redistribution_l1_j_m2+record['projection_energy_redistribution_l1_j_m2'])
        return result, advanced, ledger
