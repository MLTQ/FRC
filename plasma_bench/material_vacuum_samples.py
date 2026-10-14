"""Physical-point field reconstruction on the actual deformed magnetic exterior."""
import numpy as np
from .compatible_polynomials import bases
from .open_boundary import multipoles


def physical_samples(vacuum,coordinates,boundary,r,z,*,geometry=None):
    """Return actual finite-space psi,Br,Bz, without quadrature-cloud interpolation."""
    r,z=np.broadcast_arrays(np.asarray(r,float),np.asarray(z,float));shape=r.shape;r=r.ravel();z=z.ravel()
    if not np.isfinite(r+z).all() or np.any(r<0):raise ValueError('Finite physical exterior points with r>=0 required')
    g=vacuum.geometry(coordinates) if geometry is None else geometry;c=g['extension']@boundary;mode=c[vacuum.mode_start:]
    psi=np.empty_like(r);br=np.empty_like(r);bz=np.empty_like(r)
    outside=(r>=vacuum.outer.radius)|(abs(z)>=vacuum.outer.half_length)
    if outside.any():
        p,pr,pz=multipoles(r[outside],z[outside],vacuum.scale,vacuum.modes);rr=r[outside];zz=z[outside]
        psi[outside]=p@mode;br[outside]=np.divide(-pz@mode,rr,out=np.zeros_like(rr),where=rr>0)
        ell=np.arange(1,vacuum.modes+1);distance=np.hypot(rr,zz)
        axis=((vacuum.scale/distance[:,None])**ell/distance[:,None]**2*np.sqrt(ell*(ell+1))*np.sign(zz[:,None])**(ell-1))@mode
        bz[outside]=np.divide(pr@mode,rr,out=axis,where=rr>0)
    displacement=np.asarray(coordinates)-vacuum.reference;patches=[]
    for patch in vacuum.patches:
        rn,zn=patch['r'],patch['z'];rr,zz=np.meshgrid(rn,zn,indexing='ij')
        nodes=np.stack((rr,zz),axis=-1)+(vacuum._geometry_lift(rn,zn)@displacement).reshape(len(rn),len(zn),2)
        corners=np.stack((nodes[:-1,:-1],nodes[1:,:-1],nodes[1:,1:],nodes[:-1,1:]),axis=-2).reshape(-1,4,2)
        patches.append((patch,nodes,corners,corners.min(axis=1),corners.max(axis=1)))
    for index in np.flatnonzero(~outside):
        point=np.array([r[index],z[index]]);found=False
        for patch,nodes,corners,lower,upper in patches:
            candidates=np.flatnonzero(np.all(point>=lower-1e-13,axis=1)&np.all(point<=upper+1e-13,axis=1))
            for cell in candidates:
                quad=corners[cell];uv=np.full(2,.5)
                for _ in range(12):
                    u,v=uv;weights=np.array([(1-u)*(1-v),u*(1-v),u*v,(1-u)*v])
                    deriv=np.array([[-(1-v),-(1-u)],[1-v,-u],[v,u],[-v,1-u]])
                    jac=quad.T@deriv;error=weights@quad-point
                    if np.linalg.norm(error)<1e-13:break
                    uv-=np.linalg.solve(jac,error)
                if np.any(uv< -1e-10) or np.any(uv>1+1e-10) or np.linalg.norm(error)>1e-11:continue
                i,j=divmod(cell,len(patch['z'])-1);rn,zn=patch['r'],patch['z']
                R=rn[i]+uv[0]*(rn[i+1]-rn[i]);Z=zn[j]+uv[1]*(zn[j+1]-zn[j])
                N,D,_,_=bases(.5*rn**2,np.array([.5*R*R]),vacuum.degree)
                Nz,Dz,_,_=bases(zn,np.array([Z]),vacuum.degree)
                values=np.asarray(patch['lift']@c).reshape(len(rn),len(zn))
                cp,cr,cz,cbz=vacuum._correction(patch,np.array([R]),np.array([Z]),grid=False)
                value=float((N@values@Nz.T)[0,0]+(cp@mode)[0])
                pr=float(R*(D@values@Nz.T)[0,0]+(cr@mode)[0]);pz=float((N@values@Dz.T)[0,0]+(cz@mode)[0])
                F=jac@np.diag([1/(rn[i+1]-rn[i]),1/(zn[j+1]-zn[j])]);gradient=np.linalg.solve(F.T,np.array([pr,pz]))
                psi[index]=value
                if r[index]>0:br[index]=-gradient[1]/r[index];bz[index]=gradient[0]/r[index]
                else:br[index]=0.;bz[index]=float((D@values@Nz.T)[0,0]+(cbz@mode)[0])/F[0,0]**2
                found=True;break
            if found:break
        if not found:raise ValueError('Point lies inside the moving plasma or could not be inverted in the exterior mesh')
    return psi.reshape(shape),br.reshape(shape),bz.reshape(shape)
