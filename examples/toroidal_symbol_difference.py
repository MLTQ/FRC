"""Independent physical curl/primitive-equation check of a local toroidal symbol."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0, atomic_mass


def direct_matrix(ref,r0,z0,angle,scale,*,ion_pressure=16.,ion_mass=39.948*atomic_mass,wavenumber=1000.,spacing=2e-6):
    """Complex-step perturbations; fourth-order nested spatial differences."""
    kr,kz=np.cos(angle),np.sin(angle);epsilon=1e-20
    columns=[]
    for column in range(5):
        def fields(r,z):
            f=ref.fields(r,z)
            out=np.array([f['density'],0.,0.,ion_pressure,f['pressure'],r*f['g']],complex)
            value=1j*epsilon*scale[column]*np.sin(wavenumber*(kr*(r-r0)+kz*(z-z0)))
            if column==1:
                out[1]+=kr*value;out[2]+=kz*value
            else:
                out[(0,1,3,4,5)[column]]+=value
            return out
        def gradient(fn,r,z,axis):
            v=[fn(r+i*spacing if axis==0 else r,z+i*spacing if axis==1 else z) for i in (-2,-1,1,2)]
            return (v[0]-8*v[1]+8*v[2]-v[3])/(12*spacing)
        def electric(r,z):
            q=fields(r,z);dr,dz=gradient(fields,r,z,0),gradient(fields,r,z,1)
            n,vr,vz,pi,pe,b=q;jr=-dz[5]/MU0;jz=(dr[5]+b/r)/MU0
            return np.array([vz*b-jz*b/(n*E)-dr[4]/(n*E),-vr*b+jr*b/(n*E)-dz[4]/(n*E)])
        q=fields(r0,z0);dr,dz=gradient(fields,r0,z0,0),gradient(fields,r0,z0,1)
        n,vr,vz,pi,pe,b=q;div=dr[1]+vr/r0+dz[2]
        jr=-dz[5]/MU0;jz=(dr[5]+b/r0)/MU0
        dn=-vr*dr[0]-vz*dz[0]-n*div
        dvr=-vr*dr[1]-vz*dz[1]-(dr[3]+dr[4]+jz*b)/(ion_mass*n)
        dvz=-vr*dr[2]-vz*dz[2]-(dz[3]+dz[4]-jr*b)/(ion_mass*n)
        dpi=-vr*dr[3]-vz*dz[3]-ref.gamma*pi*div
        dpe=-vr*dr[4]-vz*dz[4]-ref.gamma*pe*div+(jr*(dr[4]-ref.gamma*pe*dr[0]/n)+jz*(dz[4]-ref.gamma*pe*dz[0]/n))/(n*E)
        db=gradient(electric,r0,z0,0)[1]-gradient(electric,r0,z0,1)[0]
        derivative=np.array([dn,kr*dvr+kz*dvz,dpi,dpe,db])
        columns.append(derivative.imag/(epsilon*wavenumber*scale))
    return np.array(columns).T
