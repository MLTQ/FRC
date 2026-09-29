"""Step-doubled resistive accuracy control without energy or state corrections."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0


@dataclass
class DiffusionResult:
    material_a: np.ndarray
    cell_heat: np.ndarray
    source_work: float
    accepted_substeps: int
    rejected_trials: int
    maximum_accepted_error: float
    maximum_balance_residual: float
    maximum_linear_residual: float


def integrate_diffusion(operator,material_a,duration,*,bias_start=0.,bias_end=None,
                        tolerance=1e-4,max_trials=10000):
    """Pure frozen-geometry integration; bias interpolates linearly over this span."""
    if not np.isfinite(duration) or duration <= 0 or not np.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError("Positive finite duration and diffusion tolerance in (0,1) required")
    bias_end=bias_start if bias_end is None else bias_end
    if not np.isfinite(bias_start) or not np.isfinite(bias_end):
        raise ValueError("Finite prescribed bias endpoints required")
    if not isinstance(max_trials,int) or max_trials < 1:
        raise ValueError("Positive integer diffusion trial budget required")
    k=operator.stiffness;m=operator.mass
    volume=operator.domain.volumes(operator.elements)
    a=np.asarray(material_a,float).copy()
    def norm(matrix,value):
        return float(np.sqrt(max(float(value@(matrix@value)),0.)))
    initial=operator.induced(a,bias_start)
    field_anchor=norm(k,initial)
    current_anchor=norm(m,operator.current_over_radius(a,bias_start))
    energy_scale=float(np.pi/MU0*(initial@(k@initial)))
    heat_anchor=energy_scale/np.sqrt(volume.sum())
    heat=np.zeros_like(volume);work=0.;elapsed=0.;h=duration
    accepted=rejected=0;maximum_error=balance=linear=0.
    tiny=np.finfo(float).tiny
    def bias(time):
        return float(bias_start+(bias_end-bias_start)*time/duration)
    for _ in range(max_trials):
        remaining=duration-elapsed
        if remaining <= 8*np.finfo(float).eps*duration:
            return DiffusionResult(a,heat,work,accepted,rejected,maximum_error,balance,linear)
        h=min(h,remaining)
        if h <= np.finfo(float).eps*max(duration,elapsed):
            raise RuntimeError("Resistive accuracy timestep underflow; accepted plasma state retained")
        b0,bhalf,b1=bias(elapsed),bias(elapsed+.5*h),bias(elapsed+h)
        coarse=operator.step(a,h,bias_start=b0,bias_end=b1)
        first=operator.step(a,.5*h,bias_start=b0,bias_end=bhalf)
        fine=operator.step(first.material_a,.5*h,bias_start=bhalf,bias_end=b1)
        qcoarse=operator.current_over_radius(coarse.material_a,b1)
        qfine=operator.current_over_radius(fine.material_a,b1)
        qheat=first.cell_heat+fine.cell_heat
        field_scale=max(field_anchor,norm(k,fine.induced),norm(k,coarse.induced),tiny)
        current_scale=max(current_anchor,norm(m,qfine),norm(m,qcoarse),tiny)
        def heat_norm(value):
            return float(np.sqrt(np.sum(value*value/volume)))
        heat_scale=max(heat_anchor,heat_norm(qheat),heat_norm(coarse.cell_heat),tiny)
        error=max(norm(k,fine.induced-coarse.induced)/field_scale,
                  norm(m,qfine-qcoarse)/current_scale,
                  heat_norm(qheat-coarse.cell_heat)/heat_scale)
        if not np.isfinite(error):
            raise ValueError("Nonfinite diffusion error estimate; accepted plasma state retained")
        if error <= tolerance:
            a=fine.material_a;heat+=qheat;work+=first.source_work+fine.source_work
            elapsed+=h;accepted+=2;maximum_error=max(maximum_error,error)
            balance=max(balance,abs(first.balance_residual),abs(fine.balance_residual))
            linear=max(linear,first.linear_relative_residual,fine.linear_relative_residual)
            field_anchor=max(field_anchor,field_scale);current_anchor=max(current_anchor,current_scale)
            heat_anchor=max(heat_anchor,heat_scale)
            factor=2. if error == 0 else min(2.,max(.2,.9*np.sqrt(tolerance/error)))
            h*=factor
        else:
            rejected+=1;h*=max(.1,min(.8,.9*np.sqrt(tolerance/error)))
    if duration-elapsed <= 8*np.finfo(float).eps*duration:
        return DiffusionResult(a,heat,work,accepted,rejected,maximum_error,balance,linear)
    raise RuntimeError("Resistive accuracy trial budget exceeded; accepted plasma state retained")
