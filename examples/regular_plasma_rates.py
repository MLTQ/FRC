"""Measured unforced rates of the actual fully assembled moving system."""
import time
import numpy as np
from .regular_plasma_reference import initialize
from .regular_plasma_norms import measure


def run_case(spec):
    start = time.monotonic()
    ref, d, u, system, field = initialize(spec)
    rate = system.derivative(u, field)
    result = dict(spec=spec, status='complete', residual=rate.residual, diagnostics=rate.diagnostics,
                  independent=measure(spec, ref, d, u, system, field, rate, 12),
                  elapsed_seconds=time.monotonic()-start)
    arrays = dict(state=u, magnetic=field.packed_increment(system.geometry), fluid_rate=rate.fluid,
                  magnetic_rate=rate.magnetic, current_rate=rate.current, electric=rate.electric)
    return result, arrays
