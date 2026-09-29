"""Independent physical interpolation, current energy and exact-wave norms."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0


def scalar(value, x, z):
    lower = (1-z)*value[:-1]+z*np.roll(value[:-1], -1, axis=1)
    upper = (1-z)*value[1:]+z*np.roll(value[1:], -1, axis=1)
    return (1-x)*lower+x*upper


def vectors(mesh, psi, g, coefficients, r, z):
    """Physical B and electric-space vector on a common point per cell."""
    x = (r*r-mesh.rf[:-1,None]**2)/(2*mesh.radial_measure[:,None])
    br = -((1-x)*np.diff(psi[:-1],axis=1)+x*np.diff(psi[1:],axis=1))/(r*mesh.dz)
    bz = ((1-z)*np.diff(psi[:,:-1],axis=0)+z*np.diff(psi[:,1:],axis=0))/mesh.radial_measure[:,None]
    b = np.stack(np.broadcast_arrays(br,r*g,bz),axis=-1)
    count = mesh.nr*mesh.nz; theta_count=(mesh.nr-1)*mesh.nz
    cr=coefficients[:count].reshape(mesh.nr,mesh.nz)
    ct=np.pad(coefficients[count:count+theta_count].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
    cz=coefficients[count+theta_count:].reshape(mesh.nr+1,mesh.nz)
    vr=r*((1-z)*cr+z*np.roll(cr,-1,axis=1))
    vt=scalar(mesh.rf[:,None]*ct,x,z)/r
    vz=(1-x)*cz[:-1]+x*cz[1:]
    return b,np.stack(np.broadcast_arrays(vr,vt,vz),axis=-1),x


def physical_energy(mesh, psi, g, current, density, ue, ui, *, reduced_mass,
                    density_scale=1e20, gamma_e=5/3, gamma_i=5/3, order=10):
    """Cell averages directly integrated in physical r,z; no production maps."""
    nodes,weights=np.polynomial.legendre.leggauss(order)
    result={k:np.zeros_like(g) for k in ('magnetic','current_kinetic','electron','ion')}
    pe=(density/density_scale)**(gamma_e-1)*ue;pi=(density/density_scale)**(gamma_i-1)*ui
    for a,wa in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*a
        for b,wb in zip(nodes,weights):
            z=.5*(1+b);mag,j,x=vectors(mesh,psi,g,current,r,z)
            weight=.25*wa*wb*r*mesh.dr/mesh.radial_measure[:,None]
            n=scalar(density,x,z)
            result['magnetic']+=weight*np.sum(mag*mag,axis=-1)/(2*MU0)
            result['current_kinetic']+=weight*reduced_mass*np.sum(j*j,axis=-1)/(2*n*E*E)
            result['electron']+=weight*scalar(pe,x,z)/(gamma_e-1)
            result['ion']+=weight*scalar(pi,x,z)/(gamma_i-1)
    return result


def wave_norms(mesh, ref, time, psi, g, current, *, order=10):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    error={k:0. for k in ('magnetic','current')}; norm=error.copy()
    axis_error=axis_norm=kinetic=exact_kinetic=0.
    for a,wa in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*a
        for b,wb in zip(nodes,weights):
            z=mesh.z[None]+.5*mesh.dz*b
            magnetic,j,_=vectors(mesh,psi,g,current,r,.5*(1+b))
            exact=ref.fields(r,z,time)
            weight=.5*np.pi*wa*wb*r*mesh.dr*mesh.dz
            for key,value in (('magnetic',magnetic),('current',j)):
                target=exact[key];scale=exact['magnetic_perturbation'] if key=='magnetic' else target
                defect=np.sum((value-target)**2,axis=-1);squared=np.sum(scale*scale,axis=-1)
                error[key]+=np.sum(weight*defect);norm[key]+=np.sum(weight*squared)
                if key=='current':
                    axis_error+=np.sum(weight[0]*defect[0]);axis_norm+=np.sum(weight[0]*squared[0])
            kinetic+=np.sum(weight*ref.masses[2]*np.sum(j*j,axis=-1)/(2*ref.number_density*E*E))
            exact_kinetic+=np.sum(weight*ref.masses[2]*np.sum(exact['current']**2,axis=-1)/(2*ref.number_density*E*E))
    return {k+'_relative_L2':float(np.sqrt(error[k]/norm[k])) for k in error} | dict(
        magnetic_rms_error_t=float(np.sqrt(error['magnetic']/np.sum(mesh.volumes))),
        current_rms_error_a_m2=float(np.sqrt(error['current']/np.sum(mesh.volumes))),
        axis_current_relative_L2=float(np.sqrt(axis_error/axis_norm)),
        physical_current_energy_j=float(kinetic),exact_current_energy_j=float(exact_kinetic))


def exact_flux(ref,r,z,time):
    f=ref.fields(r,z,time);alpha,beta,reduced=ref.masses
    kinetic=.5*ref.number_density*reduced*np.sum(f['drift']**2,axis=-1)
    enthalpy=2.5*ref.number_density*E*(alpha*ref.ion_ev-beta*ref.electron_ev)
    return np.cross(f['electric'],f['magnetic'])/MU0+f['drift']*(enthalpy-(beta-alpha)*kinetic)[...,None]


def wave_face_power(mesh,ref,time,order=12):
    """Independent continuum integrals on every radial and axial cell face."""
    nodes,weights=np.polynomial.legendre.leggauss(order)
    radial=np.zeros((mesh.nr+1,mesh.nz));axial=np.zeros((mesh.nr,mesh.nz+1))
    for a,wa in zip(nodes,weights):
        radial+=np.pi*mesh.rf[:,None]*mesh.dz*wa*exact_flux(ref,mesh.rf[:,None],mesh.z[None]+.5*mesh.dz*a,time)[...,0]
        r=mesh.r[:,None]+.5*mesh.dr*a
        axial+=np.pi*r*mesh.dr*wa*exact_flux(ref,r,mesh.zf[None],time)[...,2]
    return radial,axial
