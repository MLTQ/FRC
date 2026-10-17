"""Explicit high-precision geometry; physical coefficients remain float64."""
from decimal import Decimal, localcontext
from functools import wraps
import os
import numpy as np
from .precise_search import searchsorted


ARITHMETIC = os.environ.get('FRC_GEOMETRY_BACKEND', 'decimal')
if ARITHMETIC not in ('decimal', 'mpfr'):
    raise ValueError('FRC_GEOMETRY_BACKEND must be decimal or mpfr')
BITS = 192 if ARITHMETIC == 'mpfr' else None
DIGITS = 57 if ARITHMETIC == 'mpfr' else 50
if ARITHMETIC == 'mpfr':
    try:
        import gmpy2
    except ImportError as error:
        raise ImportError('MPFR geometry requires the plasma-bench[precision] extra') from error


def geometry_precision(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        if ARITHMETIC == 'mpfr':
            with gmpy2.context(precision=BITS, trap_divzero=True, trap_invalid=True):
                return function(*args, **kwargs)
        with localcontext() as context:
            context.prec = DIGITS
            return function(*args, **kwargs)
    return wrapped


def decimal(value):
    if ARITHMETIC == 'mpfr':
        if isinstance(value, gmpy2.mpfr):
            return value
        return gmpy2.mpfr(value if isinstance(value, str) else float(value), BITS)
    if isinstance(value, Decimal):
        return value
    return Decimal(value) if isinstance(value, str) else Decimal.from_float(float(value))


def coordinates(values):
    a = np.asarray(values)
    return np.array([decimal(x) for x in a.ravel()], object).reshape(a.shape)


@geometry_precision
def widths(values, axis=-1):
    exact = np.diff(values, axis=axis)
    result = np.asarray(exact, float)
    if np.any((exact > 0) & (result <= 0)):
        raise ValueError('Positive geometric width is outside float coefficient range')
    return result


@geometry_precision
def fraction(points, left, right):
    return np.asarray((points-left)/(right-left), float)


@geometry_precision
def interpolate(points, edges, values):
    x = coordinates(points)
    j = np.clip(searchsorted(edges, x, side='right')-1, 0, len(edges)-2)
    a = fraction(x, edges[j], edges[j+1])
    return (1-a)*values[j]+a*values[j+1]


@geometry_precision
def quadrature(edges, count):
    nodes, weights = np.polynomial.legendre.leggauss(count)
    local = (coordinates(nodes)+decimal(1))/decimal(2)
    points = edges[:-1, None]+np.diff(edges)[:, None]*local
    return points, widths(edges)[:, None]*weights/2
