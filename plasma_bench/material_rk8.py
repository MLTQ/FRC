"""Fixed DOP853 central advance of the complete material plasma state."""
import numpy as np
from scipy.integrate import DOP853


def dop853_step(model, state, dt, time=0.):
    """Use the twelve-stage eighth-order tableau without adaptive control."""
    if not np.isfinite([dt, time]).all() or dt <= 0:
        raise ValueError('Positive finite step and finite time required')
    stages = []
    for row, abscissa in zip(DOP853.A, DOP853.C):
        stage = state
        for coefficient, rate in zip(row, stages):
            if coefficient:
                stage = stage.added(rate, dt*coefficient)
        stages.append(model.derivative(stage, time=time+dt*abscissa))
    result = state
    for coefficient, rate in zip(DOP853.B, stages):
        if coefficient:
            result = result.added(rate, dt*coefficient)
    model.response(result, rate=False, time=time+dt)
    return result
