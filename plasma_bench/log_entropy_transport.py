"""Atomic RK4 verification transport in an explicitly prescribed static current."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E


@dataclass(frozen=True)
class EntropyTransportStep:
    entropy: np.ndarray
    max_scalar_residual: float
    iterations: int


class LogEntropyTransport:
    """Thermal transport fixture only: current/density are externally held fixed."""
    def __init__(self,thermal,current):
        current=np.asarray(current,float)
        if current.shape!=(thermal.model.mass.shape[0],) or not np.isfinite(current).all():
            raise ValueError('Finite fixed compatible current required')
        self.thermal=thermal;self.scalar=thermal.scalar
        self.current=current.copy();self.current.setflags(write=False)
        self.samples=thermal.geometry.sample_electric(current);self.samples.setflags(write=False)

    def rate(self,entropy):
        self.thermal.fields(entropy)
        physical=np.sum(self.samples*self.scalar.gradient(entropy),axis=-1)/E
        return self.scalar.solve(self.scalar.adjoint(physical))

    def step(self,entropy,dt):
        s=self.scalar.validate(entropy)
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite time step required')
        k1=self.rate(s);k2=self.rate(s+.5*dt*k1.values)
        k3=self.rate(s+.5*dt*k2.values);k4=self.rate(s+dt*k3.values)
        result=s+dt*(k1.values+2*k2.values+2*k3.values+k4.values)/6
        self.thermal.fields(result)
        stages=(k1,k2,k3,k4)
        return EntropyTransportStep(result,max(k.residual for k in stages),sum(k.iterations for k in stages))
