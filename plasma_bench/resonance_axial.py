"""Resonance transport with exact averaging over a uniform axial source."""
import numpy as np
from scipy.constants import c
from scipy.special import voigt_profile
from .resonance_transport import isotropic, redistribute
from .resonance_polarization import dipole_emission, observer_stokes


def resonant_decorrelation_rate(model, line, density_m3):
    """ORNL-TM-4368 p35 resonant broadening model, citing Berman & Lamb 1969."""
    return (0.904*model.A[line]*model.weight[model.upper[line]]
            /model.weight[model.lower[line]]*density_m3*model.wavelength[line]**3/(6*np.pi**2))


def axial_transport(model, line, density_m3, temperature_j, radius_m, length_m,
                    beam_radius_m, aperture_radius_m, slit_width_m, edges_s,
                    packets, seed, *, decorrelation_rate_s=0., progress=None,
                    scattering='isotropic', detector_period_s=None,
                    spectral_thresholds_relative=()):
    """Each random history represents all initial axial source positions.

    The surviving source interval is intersected with the finite end boundaries
    at every flight. Lost source intervals leave at their actual affine crossing
    times. No end loss is ignored and no surviving weight is renormalized.
    """
    if (not np.isfinite([density_m3, temperature_j, radius_m, length_m,
                        beam_radius_m, aperture_radius_m, slit_width_m, decorrelation_rate_s]).all()
            or density_m3 < 0 or decorrelation_rate_s < 0 or temperature_j <= 0
            or not 0 <= beam_radius_m < radius_m or length_m <= 0
            or not 0 < aperture_radius_m < radius_m or slit_width_m <= 0 or packets < 1):
        raise ValueError('Physical finite cylinder, beam, aperture and nonnegative rates required')
    if model.lower[line] != 0 or model.decay_rate[model.upper[line]] != model.A[line]:
        raise ValueError('Single-branch ground resonance required')
    if scattering not in ('isotropic', 'dipole'):
        raise ValueError('Explicit isotropic or dipole scattering required')
    if scattering == 'dipole' and (model.weight[0] != 1 or model.weight[model.upper[line]] != 3):
        raise ValueError('Dipole mode requires an isolated J=0-1-0 transition')
    if detector_period_s is not None and (not np.isfinite(detector_period_s) or detector_period_s <= 0):
        raise ValueError('Positive finite pulse period required')
    cuts = np.asarray(spectral_thresholds_relative, float)
    if cuts.ndim != 1 or not np.isfinite(cuts).all() or np.any(cuts <= 0):
        raise ValueError('Positive finite relative spectral thresholds required')
    edges = np.asarray(edges_s, float)
    if edges[0] != 0 or np.any(np.diff(edges) <= 0):
        raise ValueError('Increasing detector time edges beginning at zero required')
    rng = np.random.default_rng(seed); rate = model.A[line]
    sigma, gamma = model.widths(temperature_j)
    damping = (gamma[line]+decorrelation_rate_s/(4*np.pi))/sigma[line]
    strength = (density_m3*rate*model.wavelength[line]**2/(8*np.pi)
                *model.weight[model.upper[line]]/model.weight[0]/sigma[line])
    radial = beam_radius_m*np.sqrt(rng.random(packets)); phi = rng.uniform(0, 2*np.pi, packets)
    position = np.column_stack((radial*np.cos(phi), radial*np.sin(phi), np.zeros(packets)))
    initial_xy = position[:, :2].copy()
    low = np.full(packets, -length_m/2); high = -low
    clock = np.zeros(packets)
    if scattering == 'dipole':
        dipoles = isotropic(rng, packets)
        direction, polarization = dipole_emission(rng, dipoles)
    else:
        direction = isotropic(rng, packets)
    velocity = rng.normal(size=(packets, 3)); rest = damping*rng.standard_cauchy(packets)
    offset = rest+np.sum(velocity*direction, axis=1); observer_offset = rest+velocity[:, 0]
    histogram = np.zeros(len(edges)-1); integral = np.zeros(packets); outside = 0.
    stokes_histogram = np.zeros((len(histogram), 3)); stokes_integral = np.zeros((packets, 3))
    spectral_histogram = np.zeros((len(histogram), len(cuts))); spectral_integral = np.zeros(len(cuts))
    polarization_traces = []
    exits = []; traces = []; escaped_weight = np.zeros(packets); mean_exit = np.zeros(packets)
    absorption_weight = np.zeros(packets); excited_time = np.zeros(packets)
    collision_weight = np.zeros(packets); emission_weight = np.zeros(packets)
    maximum_offset = float(abs(offset).max()); maximum_view_offset = float(abs(observer_offset).max())
    alive = np.arange(packets)

    def exit_interval(ids, lo, hi, first, last, side):
        weight = (hi-lo)/length_m
        escaped_weight[ids] += weight
        mean_exit[ids] += weight*(first+last)/2
        exits.append(np.column_stack((ids, lo, hi, first, last, np.full(len(ids), side))))

    for generation in range(1_000_000):
        p = position[alive]; u = direction[alive]
        weight = (high[alive]-low[alive])/length_m
        hold = rng.exponential(1/rate, len(alive)); emission = clock[alive]+hold
        excited_time[alive] += weight*hold; emission_weight[alive] += weight
        if decorrelation_rate_s:
            reset = rng.random(len(alive)) < -np.expm1(-decorrelation_rate_s*hold)
            ids = alive[reset]
            collision_weight[ids] += weight[reset]
            if scattering == 'dipole':
                dipoles[ids] = isotropic(rng, len(ids))
                direction[ids], polarization[ids] = dipole_emission(rng, dipoles[ids])
                u[reset] = direction[ids]
            v = rng.normal(size=(len(ids), 3)); rest = damping*rng.standard_cauchy(len(ids))
            offset[ids] = rest+np.sum(v*u[reset], axis=1)
            observer_offset[ids] = rest+v[:, 0]
        maximum_offset = max(maximum_offset, float(abs(offset[alive]).max()))
        maximum_view_offset = max(maximum_view_offset, float(abs(observer_offset[alive]).max()))

        # Exact initial-source measure visible through the slit at this event.
        visible = abs(p[:, 1]) < aperture_radius_m
        cap = np.zeros(len(alive))
        cap[visible] = np.minimum(slit_width_m/2,
            np.sqrt(aperture_radius_m**2-p[visible, 1]**2))
        viewed_low = np.maximum(low[alive], -p[:, 2]-cap)
        viewed_high = np.minimum(high[alive], -p[:, 2]+cap)
        viewed_length = np.maximum(viewed_high-viewed_low, 0.)
        distance_view = np.sqrt(radius_m**2-p[:, 1]**2)-p[:, 0]
        transmission = viewed_length/length_m*np.exp(
            -strength*voigt_profile(observer_offset[alive], 1., damping)*distance_view)
        angular = observer_stokes(dipoles[alive]) if scattering == 'dipole' else np.tile([1/(4*np.pi), 0., 0.], (len(alive), 1))
        stokes = transmission[:, None]*angular; score = stokes[:, 0]
        arrival = emission+distance_view/c
        detector_time = arrival if detector_period_s is None else np.remainder(arrival, detector_period_s)
        index = np.searchsorted(edges, detector_time, side='right')-1
        in_bins = (index >= 0) & (index < len(histogram))
        histogram += np.bincount(index[in_bins], weights=score[in_bins], minlength=len(histogram))
        integral[alive] += score; outside += float(score[~in_bins].sum())
        stokes_integral[alive] += stokes
        for component in range(3):
            stokes_histogram[:, component] += np.bincount(index[in_bins], weights=stokes[in_bins, component], minlength=len(histogram))
        relative_frequency = abs(observer_offset[alive])*sigma[line]/model.frequency[line]
        for j, cut in enumerate(cuts):
            wing_score = score*(relative_frequency > cut)
            spectral_integral[j] += wing_score.sum()
            spectral_histogram[:, j] += np.bincount(index[in_bins], weights=wing_score[in_bins], minlength=len(histogram))

        a = np.sum(u[:, :2]**2, axis=1); b = np.sum(p[:, :2]*u[:, :2], axis=1)
        inside = radius_m**2-np.sum(p[:, :2]**2, axis=1)
        root = np.sqrt(b*b+a*inside)
        side_distance = np.full(len(alive), np.inf)
        np.divide(inside, root+b, out=side_distance, where=(b >= 0) & (root+b > 0))
        np.divide(root-b, a, out=side_distance, where=(b < 0) & (a > 0))
        opacity = strength*voigt_profile(offset[alive], 1., damping)
        free_path = np.full(len(alive), np.inf)
        np.divide(rng.exponential(size=len(alive)), opacity, out=free_path, where=opacity > 0)
        distance = np.minimum(free_path, side_distance)
        if not np.isfinite(distance).all() or np.any(distance < 0):
            raise ValueError('Nonfinite/negative radial flight proposal')
        end = emission+distance/c; moved = p+distance[:, None]*u
        # These min/max operations intersect geometric source intervals; they do
        # not clip plasma state, frequency, photon weight or an escaped packet.
        new_low = np.minimum(np.maximum(low[alive], -length_m/2-moved[:, 2]), high[alive])
        new_high = np.maximum(np.minimum(high[alive], length_m/2-moved[:, 2]), low[alive])
        lost_low = np.where(u[:, 2] > 0, new_high, low[alive])
        lost_high = np.where(u[:, 2] > 0, high[alive], new_low)
        lost = lost_high > lost_low
        if np.any(lost):
            numerator = np.sign(u[lost, 2])*length_m/2-p[lost, 2]
            t0 = emission[lost]+(numerator-lost_low[lost])/(c*u[lost, 2])
            t1 = emission[lost]+(numerator-lost_high[lost])/(c*u[lost, 2])
            exit_interval(alive[lost], lost_low[lost], lost_high[lost],
                          np.minimum(t0,t1), np.maximum(t0,t1), 0)
        survives = new_high > new_low
        radial_exit = survives & (free_path >= side_distance)
        if np.any(radial_exit):
            exit_interval(alive[radial_exit], new_low[radial_exit], new_high[radial_exit],
                          end[radial_exit], end[radial_exit], 1)
        low[alive] = new_low; high[alive] = new_high
        position[alive] = moved; clock[alive] = end
        for local in np.flatnonzero(alive < 4):
            traces.append([int(alive[local]), generation, float(emission[local]), float(end[local]),
                           *moved[local], float(new_low[local]), float(new_high[local]), float(offset[alive[local]])])
            if scattering == 'dipole':
                k = alive[local]
                polarization_traces.append([int(k), generation, *dipoles[k], *direction[k], *polarization[k]])
        alive = alive[survives & ~radial_exit]
        if progress is not None and generation % 1000 == 0:
            progress(generation, len(alive))
        if not len(alive):
            break
        absorption_weight[alive] += (high[alive]-low[alive])/length_m
        if scattering == 'dipole':
            dipoles[alive] = polarization[alive]
            outgoing, polarization[alive] = dipole_emission(rng, dipoles[alive])
        else:
            outgoing = isotropic(rng, len(alive))
        offset[alive], observer_offset[alive] = redistribute(
            rng, offset[alive], direction[alive], outgoing, damping, observer=True)
        direction[alive] = outgoing
    else:
        raise RuntimeError('Live axial source remains; generation limit is not escape')
    return dict(detector_time_edges_s=edges, detector_probability_per_sr=histogram/packets,
        scattering=scattering, detector_period_s=np.nan if detector_period_s is None else detector_period_s,
        detector_stokes_per_sr=stokes_histogram/packets, detector_packet_stokes_integral_per_sr=stokes_integral,
        spectral_thresholds_relative=cuts, detector_spectral_tail_probability_per_sr=spectral_histogram/packets,
        detector_spectral_tail_integral_per_sr=spectral_integral/packets,
        example_polarization_events=np.asarray(polarization_traces),
        detector_packet_integral_per_sr=integral, detector_outside_time_probability_per_sr=outside/packets,
        escaped_source_weight=escaped_weight, mean_exit_time_per_history_s=mean_exit,
        absorption_source_weight=absorption_weight, excited_residence_per_history_s=excited_time,
        collision_emission_source_weight=collision_weight, emission_source_weight=emission_weight,
        exit_intervals=np.concatenate(exits), initial_xy_m=initial_xy, example_events=np.asarray(traces),
        maximum_relative_detuning=maximum_offset*sigma[line]/model.frequency[line],
        maximum_relative_observer_detuning=maximum_view_offset*sigma[line]/model.frequency[line],
        collision_probability=decorrelation_rate_s/(rate+decorrelation_rate_s),
        decorrelation_rate_s=decorrelation_rate_s, voigt_damping_sigma=damping)
