"""Analog resonance-photon transport through a finite uniform neutral cylinder."""
import numpy as np
from scipy.constants import c
from scipy.special import ndtr, ndtri, voigt_profile


def isotropic(rng, count):
    mu = rng.uniform(-1., 1., count)
    azimuth = rng.uniform(0., 2*np.pi, count)
    radial = np.sqrt(1-mu*mu)
    return np.column_stack((radial*np.cos(azimuth), radial*np.sin(azimuth), mu))


def cylinder_distance(position, direction, radius, half_length):
    """Exact positive ray distance to the side or either end of the cylinder."""
    a = np.sum(direction[:, :2]**2, axis=1)
    b = np.sum(position[:, :2]*direction[:, :2], axis=1)
    inside = radius**2-np.sum(position[:, :2]**2, axis=1)
    root = np.sqrt(b*b+a*inside)
    side = np.full(len(a), np.inf)
    # Stable quadratic root for outgoing rays near the side wall.
    positive = b >= 0
    np.divide(inside, root+b, out=side, where=positive & (root+b > 0))
    np.divide(root-b, a, out=side, where=~positive & (a > 0))
    dz = direction[:, 2]
    end = np.full(len(a), np.inf)
    np.divide(half_length-np.sign(dz)*position[:, 2], abs(dz), out=end, where=dz != 0)
    distance = np.minimum(side, end)
    if not np.isfinite(distance).all() or np.any(distance < 0):
        raise ValueError('A photon left the cylinder before the declared boundary event')
    return distance, side <= end


def absorbing_velocity(rng, detuning, damping):
    """Exact rejection sampling of u_parallel conditional on absorption.

    u is in thermal velocity standard deviations. Density is proportional to
    exp(-u²/2) / ((x-u)²+a²), with x in Doppler standard deviations.
    """
    x = abs(np.asarray(detuning, float)); result = np.empty_like(x)
    remaining = np.arange(len(x))
    while len(remaining):
        xx = x[remaining]; sample = np.empty_like(xx); probability = np.empty_like(xx)
        core = xx <= 2.
        sample[core] = xx[core]+damping*rng.standard_cauchy(np.sum(core))
        probability[core] = np.exp(-sample[core]**2/2)
        wing = ~core
        if np.any(wing):
            xxw = xx[wing]; split = np.minimum(.8*xxw, 4.)
            normal_mass = ndtr(split)
            lorentz_bound = damping/(np.pi*((xxw-split)**2+damping**2))
            tail_mass = .5+np.arctan((xxw-split)/damping)/np.pi
            normal_bound = np.exp(-split**2/2)/np.sqrt(2*np.pi)
            weight_left = lorentz_bound*normal_mass
            weight_right = normal_bound*tail_mass
            left = rng.random(len(xxw)) < weight_left/(weight_left+weight_right)
            candidate = np.empty_like(xxw); accept = np.empty_like(xxw)
            # Left proposal: normal distribution restricted to u <= split.
            candidate[left] = ndtri((1-rng.random(np.sum(left)))*normal_mass[left])
            accept[left] = (((xxw[left]-split[left])**2+damping**2)
                            /((xxw[left]-candidate[left])**2+damping**2))
            # Right proposal: Lorentz distribution restricted to u > split.
            angle = np.pi*(.5-(1-rng.random(np.sum(~left)))*tail_mass[~left])
            candidate[~left] = xxw[~left]+damping*np.tan(angle)
            accept[~left] = np.exp((split[~left]**2-candidate[~left]**2)/2)
            sample[wing] = candidate; probability[wing] = accept
        if np.any(probability > 1+1e-13):
            raise ValueError('Invalid conditional velocity rejection envelope')
        accepted = rng.random(len(xx)) < probability
        result[remaining[accepted]] = sample[accepted]
        remaining = remaining[~accepted]
    return np.where(np.asarray(detuning) < 0, -result, result)


def redistribute(rng, detuning, incoming, outgoing, damping, *, observer=False):
    """Coherent atom-frame frequency, with conditional Maxwellian atom velocity."""
    parallel = absorbing_velocity(rng, detuning, damping)
    reference = np.zeros_like(incoming)
    reference[:, 2] = 1.
    near_pole = abs(incoming[:, 2]) > .9
    reference[near_pole] = [1., 0., 0.]
    e1 = np.cross(incoming, reference)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(incoming, e1)
    transverse = rng.normal(size=(len(detuning), 2))
    velocity = (parallel[:, None]*incoming+transverse[:, 0, None]*e1
                +transverse[:, 1, None]*e2)
    rest = detuning-parallel
    shifted = rest+np.sum(velocity*outgoing, axis=1)
    return (shifted, rest+velocity[:, 0]) if observer else shifted


def transport(model, line, density_m3, temperature_j, radius_m, length_m,
              packets, seed, sample_times_s, *, redistribution='atom_frame',
              maximum_generations=1_000_000, progress=None,
              initial_positions_m=None, side_view=None):
    """Track every weak-pulse quantum from a specified atom until wall escape.

    Each packet starts as an excited atom. The fixed neutral bath, isotropic
    scattering, exponential atom lifetimes, and narrow-line constant energy
    are explicit approximations; no photon is discarded or forced to escape.
    """
    if redistribution not in ('atom_frame', 'complete'):
        raise ValueError('Unknown redistribution prescription')
    if (not np.isfinite([density_m3, temperature_j, radius_m, length_m]).all()
            or density_m3 < 0 or min(temperature_j, radius_m, length_m) <= 0
            or packets <= 0):
        raise ValueError('Physical geometry, temperature and nonnegative density required')
    if model.lower[line] != 0 or model.decay_rate[model.upper[line]] != model.A[line]:
        raise ValueError('This transport component requires a single-branch ground resonance')
    grid = np.asarray(sample_times_s, float)
    if grid[0] != 0 or np.any(np.diff(grid) <= 0):
        raise ValueError('Strictly ordered output times beginning at zero required')
    rng = np.random.default_rng(seed)
    sigma, gamma = model.widths(temperature_j)
    damping = gamma[line]/sigma[line]
    strength = (density_m3*model.A[line]*model.wavelength[line]**2/(8*np.pi)
                *model.weight[model.upper[line]]/model.weight[0]/sigma[line])
    rate = model.A[line]
    position = (np.zeros((packets, 3)) if initial_positions_m is None
                else np.array(initial_positions_m, dtype=float, copy=True))
    if (position.shape != (packets, 3) or not np.isfinite(position).all()
            or np.any(np.sum(position[:, :2]**2, axis=1) >= radius_m**2)
            or np.any(abs(position[:, 2]) >= length_m/2)):
        raise ValueError('Every initial excitation must be strictly inside the cylinder')
    initial_position = position.copy(); clock = np.zeros(packets)
    direction = isotropic(rng, packets)
    if side_view is None:
        offset = rng.normal(size=packets)+damping*rng.standard_cauchy(packets)
    else:
        velocity = rng.normal(size=(packets, 3))
        rest = damping*rng.standard_cauchy(packets)
        offset = rest+np.sum(velocity*direction, axis=1)
        observer_offset = rest+velocity[:, 0]
    absorptions = np.zeros(packets, dtype=np.int64)
    flight_total = np.zeros(packets); exit_time = np.empty(packets)
    exit_offset = np.empty(packets); exit_side = np.zeros(packets, bool)
    first_flight_time = np.zeros(packets)
    atom_delta = np.zeros(len(grid)+1, dtype=np.int64)
    photon_delta = atom_delta.copy()
    alive = np.arange(packets); traces = []
    maximum_relative_detuning = float(abs(offset).max()*sigma[line]/model.frequency[line])
    def interval(delta, begin, end):
        start_index = np.searchsorted(grid, begin, side='left')
        end_index = np.searchsorted(grid, end, side='left')
        delta += np.bincount(start_index, minlength=len(delta))
        delta -= np.bincount(end_index, minlength=len(delta))
    for generation in range(maximum_generations):
        start = clock[alive].copy()
        emission = start+rng.exponential(1/rate, len(alive))
        interval(atom_delta, start, emission)
        if side_view is not None:
            side_view.record(alive, position[alive], emission, observer_offset[alive],
                             strength, damping)
        boundary, side = cylinder_distance(position[alive], direction[alive], radius_m, length_m/2)
        opacity = strength*voigt_profile(offset[alive], 1., damping)
        path = np.full(len(alive), np.inf)
        np.divide(rng.exponential(size=len(alive)), opacity, out=path, where=opacity > 0)
        escaped = path >= boundary
        distance = np.minimum(path, boundary)
        end = emission+distance/c
        interval(photon_delta, emission, end)
        flight_total[alive] += distance/c
        position[alive] += distance[:, None]*direction[alive]
        clock[alive] = end
        if generation == 0:
            first_flight_time[:] = distance/c
        for local in np.flatnonzero(alive < 8):
            traces.append([int(alive[local]), generation, float(emission[local]),
                           float(end[local]), *position[alive[local]].tolist(),
                           float(offset[alive[local]]), int(escaped[local])])
        ids = alive[escaped]
        exit_time[ids] = end[escaped]; exit_offset[ids] = offset[ids]
        exit_side[ids] = side[escaped]
        alive = alive[~escaped]
        if progress is not None and generation % 1000 == 0:
            progress(generation, len(alive))
        if not len(alive):
            break
        absorptions[alive] += 1
        outgoing = isotropic(rng, len(alive))
        if redistribution == 'atom_frame':
            shifted = redistribute(rng, offset[alive], direction[alive], outgoing,
                                   damping, observer=side_view is not None)
            if side_view is None:
                offset[alive] = shifted
            else:
                offset[alive], observer_offset[alive] = shifted
        else:
            if side_view is None:
                offset[alive] = rng.normal(size=len(alive))+damping*rng.standard_cauchy(len(alive))
            else:
                velocity = rng.normal(size=(len(alive), 3))
                rest = damping*rng.standard_cauchy(len(alive))
                offset[alive] = rest+np.sum(velocity*outgoing, axis=1)
                observer_offset[alive] = rest+velocity[:, 0]
        direction[alive] = outgoing
        maximum_relative_detuning = max(maximum_relative_detuning,
            float(abs(offset[alive]).max()*sigma[line]/model.frequency[line]))
    else:
        raise RuntimeError(f'{len(alive)} packets remain; generation limit is not escape')
    atoms = np.cumsum(atom_delta[:-1]); photons = np.cumsum(photon_delta[:-1])
    escaped_count = np.searchsorted(np.sort(exit_time), grid, side='right')
    return dict(time_s=grid, excited_count=atoms, photon_count=photons,
        escaped_count=escaped_count, packet_balance=atoms+photons+escaped_count-packets,
        exit_time_s=exit_time, exit_position_m=position, exit_detuning_sigma=exit_offset,
        initial_position_m=initial_position, exit_direction=direction,
        exit_side=exit_side, absorptions=absorptions, flight_time_s=flight_total,
        first_flight_time_s=first_flight_time, example_events=np.array(traces),
        maximum_relative_detuning=np.array(maximum_relative_detuning))
