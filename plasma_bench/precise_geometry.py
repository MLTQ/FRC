"""Explicit decimal geometry; physical coefficients remain float64."""
from decimal import Decimal, localcontext
from functools import wraps
import numpy as np


DIGITS = 50


def geometry_precision(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        with localcontext() as context:
            context.prec = DIGITS
            return function(*args, **kwargs)
    return wrapped


def decimal(value):
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
    j = np.clip(np.searchsorted(edges, x, side='right')-1, 0, len(edges)-2)
    a = fraction(x, edges[j], edges[j+1])
    return (1-a)*values[j]+a*values[j+1]


@geometry_precision
def quadrature(edges, count):
    nodes, weights = np.polynomial.legendre.leggauss(count)
    local = (coordinates(nodes)+Decimal(1))/Decimal(2)
    points = edges[:-1, None]+np.diff(edges)[:, None]*local
    return points, widths(edges)[:, None]*weights/2
