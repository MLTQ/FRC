"""Analytic curl-kernel separation for ill-conditioned electric backgrounds."""
import numpy as np
from .inertial_induction import InertialInduction,InductionRate
from .moving_current import MovingCurrent


def separate(geometry,electric):
    """Er/r independent of z and Ez independent of r are exact curl kernels."""
    mesh=geometry.mesh;r,t,z=geometry.split_electric(electric)
    if z.size!=(mesh.nr+1)*mesh.nz:raise ValueError('Free full radial Ez trace required')
    r=r.reshape(mesh.nr,mesh.nz);z=z.reshape(mesh.nr+1,mesh.nz)
    kr=np.broadcast_to(np.mean(r,axis=1)[:,None],r.shape)
    kz=np.broadcast_to(z[0],z.shape)
    kernel=np.concatenate((kr.ravel(),np.zeros_like(t),kz.ravel()))
    return electric-kernel,kernel


def difference_curl(geometry,electric):
    """Difference before metric division; analytic kernel vanishes bit for bit."""
    mesh=geometry.mesh;r,t,z=geometry.split_electric(electric)
    r=r.reshape(mesh.nr,mesh.nz);t=t.reshape(mesh.nr-1,mesh.nz);z=z.reshape(mesh.nr+1,mesh.nz)
    p=np.zeros((mesh.nr+1,mesh.nz));p[1:-1]=-mesh.rf[1:-1,None]*t
    g=(r-np.roll(r,-1,axis=1))/mesh.dz+np.diff(z,axis=0)/mesh.radial_measure[:,None]
    return np.r_[p.ravel(),g.ravel()]


class KernelInduction(InertialInduction):
    def solve(self,forcing):
        model=self.model;f=np.asarray(forcing,float)
        if f.shape!=(model.mass.shape[0],) or not np.isfinite(f).all():raise ValueError('Finite full weak electric forcing required')
        electric0=self.mass_inverse(f);small,kernel=separate(model.geometry,electric0)
        result=super().solve(model.mass@small)
        electric=result.electric+kernel
        norm=lambda a:float(np.linalg.norm(a))
        relative=lambda a,b:norm(a)/max(norm(b),np.finfo(float).tiny)
        # These are backward errors of the ORIGINAL full equations. The inner
        # solve separately retains every small-signal relative residual check.
        inertia=model.inertia@result.current
        ohm=relative(model.mass@electric-inertia-f,abs(model.mass@electric)+abs(inertia)+abs(f))
        curl=difference_curl(model.geometry,electric)
        bound=abs(model.geometry.faraday)@abs(electric)+abs(result.magnetic)
        faraday=relative(curl-result.magnetic,bound)
        reconstruction=relative(electric0-small-kernel,abs(electric0)+abs(small)+abs(kernel))
        components=dict(result.residual_components,full_ohm_backward=ohm,
            full_faraday_backward=faraday,kernel_reconstruction_backward=reconstruction)
        residual=max(components.values())
        if not np.isfinite(residual) or residual>self.tolerance:
            raise RuntimeError(f'Kernel-separated original equation backward errors {components}; no input modified')
        # Full-field rate-relative replay can be ill-conditioned and is reported
        # distinctly, never substituted for the inner small-signal accuracy.
        components['full_faraday_rate_relative']=relative(curl-result.magnetic,result.magnetic)
        components['kernel_to_remainder_ratio']=norm(kernel)/norm(small) if norm(small)>0 else None
        return InductionRate(result.magnetic,result.current,electric,result.iterations,residual,result.refinements,components)


class KernelMovingCurrent(MovingCurrent):
    def __init__(self,model,**solver_options):
        self.model=model;self.solver=KernelInduction(model,**solver_options)
