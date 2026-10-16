"""A narrow side-view radiance estimator for finite-cylinder photon histories."""
import numpy as np
from scipy.constants import c
from scipy.special import voigt_profile


class SideView:
    """Parallel +x viewing through a circular aperture and an axial slit.

    Every emission contributes its probability density per steradian of reaching
    the viewing aperture without absorption. The estimator does not terminate,
    redirect, split or remove any analog photon from the transport calculation.
    """
    def __init__(self, radius_m, aperture_radius_m, slit_width_m, time_edges_s,
                 source_weights):
        self.radius = radius_m
        self.aperture = aperture_radius_m
        self.half_width = slit_width_m/2
        self.edges = np.asarray(time_edges_s, dtype=float)
        self.weights = np.asarray(source_weights, dtype=float)
        if (not 0 < aperture_radius_m < radius_m or slit_width_m <= 0
                or self.edges.ndim != 1 or len(self.edges) < 2
                or np.any(np.diff(self.edges) <= 0) or self.edges[0] != 0
                or not np.isfinite(self.weights).all() or np.any(self.weights <= 0)):
            raise ValueError('Physical aperture, ordered nonnegative bins and positive source weights required')
        self.histogram = np.zeros(len(self.edges)-1)
        self.integral = np.zeros(len(self.weights))
        self.outside_time_bins = 0.
        self.maximum_relative_view_offset = 0.

    def record(self, ids, position, emission, detuning, strength, damping):
        visible = ((abs(position[:, 2]) <= self.half_width)
                   & (np.sum(position[:, 1:]**2, axis=1) <= self.aperture**2))
        if not np.any(visible):
            return
        p = position[visible]; packet_ids = ids[visible]
        distance = np.sqrt(self.radius**2-p[:, 1]**2)-p[:, 0]
        if np.any(distance < 0):
            raise ValueError('Detector emission point is outside the radial wall')
        offsets = detuning[visible]
        self.maximum_relative_view_offset = max(self.maximum_relative_view_offset,
                                                float(abs(offsets).max()))
        probability = self.weights[packet_ids]*np.exp(
            -strength*voigt_profile(offsets, 1., damping)*distance)/(4*np.pi)
        arrival = emission[visible]+distance/c
        index = np.searchsorted(self.edges, arrival, side='right')-1
        inside = (index >= 0) & (index < len(self.histogram))
        self.histogram += np.bincount(index[inside], weights=probability[inside],
                                      minlength=len(self.histogram))
        self.integral[packet_ids] += probability
        self.outside_time_bins += float(probability[~inside].sum())

    def result(self):
        return dict(detector_time_edges_s=self.edges,
                    detector_probability_per_sr=self.histogram/len(self.weights),
                    detector_packet_integral_per_sr=self.integral,
                    detector_outside_time_probability_per_sr=self.outside_time_bins/len(self.weights),
                    detector_maximum_detuning_sigma=self.maximum_relative_view_offset,
                    source_weights=self.weights)
