"""Independent finite-difference energy derivative and moving-density work audit."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .current_thermal_norms import vectors,scalar


def verify(model,psi,g,ue,ui,nt,induction,forcing):
    geo=model.geometry;packed=geo.pack(psi,g);j=model.current(packed)
    samples=geo.sample_electric(j);fraction=geo.scalar_samples(nt)/model.density
    density_work=-.5*np.sum(geo.weight*model.inertia_coefficient*np.sum(samples*samples,axis=-1)*fraction)
    magnetic=float(packed@(model.hessian@induction.magnetic))
    current=float(j@(model.inertia@induction.current))
    electric_work=float(j@forcing)
    residual=magnetic+current+electric_work
    scale=abs(magnetic)+abs(current)+abs(density_work)+abs(electric_work)
    algebraic=abs(residual)/max(scale,np.finfo(float).tiny)
    assert algebraic<1e-9
    # Evaluate the exact symmetric energy difference before quadrature, using
    # polarization/rational common denominators to avoid subtracting large totals.
    dp,dg=geo.unpack(induction.magnetic);nodes,weights=np.polynomial.legendre.leggauss(8)
    def centered_rate(h):
        total=0.;mesh=model.mesh
        for x,wx in zip(nodes,weights):
            r=mesh.r[:,None]+.5*mesh.dr*x
            for z,wz in zip(nodes,weights):
                az=.5*(1+z);b,c,ax=vectors(mesh,psi,g,j,r,az)
                bd,cd,_=vectors(mesh,dp,dg,induction.current,r,az)
                n=scalar(model.electron.density,ax,az);nd=scalar(nt,ax,az)
                numerator=2*n*np.sum(c*cd,axis=-1)-nd*(np.sum(c*c,axis=-1)+h*h*np.sum(cd*cd,axis=-1))
                difference=np.sum(b*bd,axis=-1)/MU0+model.reduced_mass/(2*E*E)*numerator/(n*n-h*h*nd*nd)
                total+=np.sum(.5*np.pi*wx*wz*r*mesh.dr*mesh.dz*difference)
        return float(total)
    # h is an audit differentiation interval, not an integration time step.
    density_scale=np.max(abs(nt)/model.electron.density)
    magnetic_scale=np.sqrt((induction.magnetic@(model.hessian@induction.magnetic))/(packed@(model.hessian@packed)))
    h=1e-3/max(density_scale,magnetic_scale,1.)
    finite=[]
    for step in (h,h/2):
        derivative=(4*centered_rate(step/2)-centered_rate(step))/3
        finite.append(dict(interval_s=step,rate_w=derivative,
            relative_balance=abs(derivative+electric_work-density_work)/max(scale,np.finfo(float).tiny)))
        assert finite[-1]['relative_balance']<2e-8
    return dict(magnetic_rate_w=magnetic,current_rate_at_fixed_density_w=current,
        changing_density_current_energy_w=float(density_work),electric_forcing_work_w=electric_work,
        mixed_equation_power_relative=algebraic,independent_physical_energy_derivative=finite)
