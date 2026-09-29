"""Stable positive logarithmic means for discrete electron thermodynamics."""
import numpy as np


def logarithmic_mean(left, right):
    left, right = np.broadcast_arrays(np.asarray(left, float), np.asarray(right, float))
    if np.any(left <= 0) or np.any(right <= 0) or not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError('Positive finite arguments required for logarithmic thermodynamic mean')
    low, high = np.minimum(left, right), np.maximum(left, right)
    average = low+.5*(high-low)
    z = .5*(high-low)/average
    small = z < 1e-3
    result = np.empty_like(average)
    z2 = z[small]**2
    result[small] = average[small]*(1-z2/3-4*z2*z2/45-44*z2*z2*z2/945)
    medium = (~small) & (z < .5)
    result[medium] = average[medium]*z[medium]/np.arctanh(z[medium])
    large = ~(small | medium)
    with np.errstate(over='ignore'):
        ratio = (high[large]-low[large])/low[large]
    log_ratio = np.log1p(ratio)
    overflow = ~np.isfinite(log_ratio)
    log_ratio[overflow] = np.log(high[large][overflow])-np.log(low[large][overflow])
    result[large] = (high[large]-low[large])/log_ratio
    return result


def entropy_mean(left, right):
    """Reciprocal logarithmic mean of reciprocals; wL*wR/L(wL,wR)."""
    low, high = np.minimum(left, right), np.maximum(left, right)
    return low*(high/logarithmic_mean(low, high))
