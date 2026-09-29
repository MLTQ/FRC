"""Independent physical energy and entropy quadrature for conservative states."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .current_thermal_norms import vectors,scalar


def audit(mesh,state,parameters,evaluation,*,order=10):
    u=state.conserved;rho=u[...,0];mass=parameters.masses
    # Independent quarter-mass overlap, including both periodic neighboring cells.
    quarter=.25*mesh.volumes*rho/mass.total;pair=quarter+np.roll(quarter,1,axis=1)
    particle=np.zeros((mesh.nr+1,mesh.nz));particle[:-1]+=pair;particle[1:]+=pair
    widths=np.diff(np.r_[0.,.25*(mesh.rf[:-1]**2+mesh.rf[1:]**2),.5*mesh.radius**2])
    dual_volume=2*np.pi*widths[:,None]*mesh.dz;nnode=particle/dual_volume
    np.testing.assert_allclose(nnode,evaluation.model.electron.density,rtol=3e-15)
    pe=(nnode/parameters.density_scale)**(parameters.gamma_e-1)*state.electron_entropy
    pi=(parameters.gamma_i-1)*evaluation.components['ion']
    r2=.5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:,None]
    vr=u[...,1]/rho;vz=u[...,2]/rho;omega=u[...,4]/(rho*r2)
    energy={k:0. for k in ('magnetic','current','bulk','electron','ion')};se=si=ne=ni=0.
    nodes,weights=np.polynomial.legendre.leggauss(order)
    for x,wx in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*x
        for z,wz in zip(nodes,weights):
            az=.5*(z+1);b,j,ax=vectors(mesh,state.psi,state.g,evaluation.current,r,az)
            n=scalar(nnode,ax,az);pressure=scalar(pe,ax,az);ion_n=rho/mass.total
            weight=.5*np.pi*wx*wz*r*mesh.dr*mesh.dz
            values=dict(magnetic=np.sum(b*b,axis=-1)/(2*MU0),current=mass.reduced*np.sum(j*j,axis=-1)/(2*n*E**2),
                bulk=.5*rho*(vr*vr+vz*vz+r*r*omega*omega),electron=pressure/(parameters.gamma_e-1),ion=pi/(parameters.gamma_i-1))
            for key,value in values.items():energy[key]+=float(np.sum(weight*value))
            se+=float(np.sum(weight*n*np.log(pressure/(n*E)*(n/parameters.density_scale)**(1-parameters.gamma_e))))
            si+=float(np.sum(weight*ion_n*np.log(pi/(ion_n*E)*(ion_n/parameters.density_scale)**(1-parameters.gamma_i))))
            ne+=float(np.sum(weight*n));ni+=float(np.sum(weight*ion_n))
    total=float(np.sum(mesh.volumes*u[...,3]));relative=abs(sum(energy.values())-total)/abs(total)
    assert relative<1e-10,relative
    assert abs(ne-ni)/ni<3e-14
    return dict(energies_j=energy,total_relative_difference=relative,electron_particles=ne,ion_particles=ni,
        electron_physical_log_entropy_particles=se,ion_physical_log_entropy_particles=si)
