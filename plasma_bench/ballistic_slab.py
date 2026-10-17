"""Independent exact Gaussian characteristics integrated over spatial cells."""
import numpy as np
from scipy.special import ndtr


def cell_moments(edges, half_width, thermal_speed, time, sigma):
    """Return cell averages of integral s**k phi(s) 1(|x-vT*s*t|<a) ds, k=0..3."""
    if time <= 0:
        raise ValueError('Characteristic endpoint requires positive time')
    travel = thermal_speed*time
    lo_x, hi_x = edges[:-1, None], edges[1:, None]
    cuts = np.sort(np.concatenate((np.full_like(lo_x, -sigma),
        np.clip((lo_x-half_width)/travel, -sigma, sigma),
        np.clip((hi_x-half_width)/travel, -sigma, sigma),
        np.clip((lo_x+half_width)/travel, -sigma, sigma),
        np.clip((hi_x+half_width)/travel, -sigma, sigma),
        np.full_like(lo_x, sigma)), axis=1), axis=1)
    lo, hi = cuts[:, :-1], cuts[:, 1:]
    def overlap(s):
        return np.maximum(0., np.minimum(hi_x, half_width+travel*s)
                          -np.maximum(lo_x, -half_width+travel*s))
    slope = np.divide(overlap(hi)-overlap(lo), hi-lo,
                      out=np.zeros_like(lo), where=hi > lo)
    intercept = overlap(lo)-slope*lo
    phi_lo, phi_hi = np.exp(-lo*lo/2)/np.sqrt(2*np.pi), np.exp(-hi*hi/2)/np.sqrt(2*np.pi)
    integrals = [np.where(lo >= 0, ndtr(-lo)-ndtr(-hi), ndtr(hi)-ndtr(lo)), phi_lo-phi_hi]
    for k in range(2, 5):
        integrals.append((k-1)*integrals[k-2]+lo**(k-1)*phi_lo-hi**(k-1)*phi_hi)
    return np.stack([np.sum(intercept*integrals[k]+slope*integrals[k+1], axis=1)/np.diff(edges)
                     for k in range(4)], axis=1)
