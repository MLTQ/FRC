"""Magnetic-only exterior around a finite cylinder, with exact outer trace matching."""
from types import SimpleNamespace
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse import csr_matrix, diags, kron
from .compatible_polynomials import bases
from .open_boundary import multipoles, exterior_energy


class FiniteCylinderVacuum:
    """Static vacuum elimination; boundary unknowns are the actual plasma flux."""
    def __init__(self, mesh, *, degree=5, quadrature_order=12, vacuum_cells=8,
                 modes=8, outer_radius=None, outer_half_length=None):
        self.mesh=mesh; self.degree=degree; self.order=quadrature_order; self.modes=modes
        if type(vacuum_cells) is not int or vacuum_cells<degree or type(modes) is not int or modes<1:
            raise ValueError('Resolved vacuum patches and positive multipole count required')
        enclosing=np.hypot(mesh.radius,mesh.half_length)
        rout=1.5*enclosing if outer_radius is None else float(outer_radius)
        zout=2*enclosing if outer_half_length is None else float(outer_half_length)
        if not np.isfinite([rout,zout]).all() or min(rout,zout)<=enclosing:
            raise ValueError('Outer matching rectangle must enclose the entire plasma sphere')
        self.outer=SimpleNamespace(radius=rout,half_length=zout); self.scale=min(rout,zout)
        rs=np.geomspace(mesh.radius,rout,vacuum_cells+1)
        zb=np.linspace(-zout,-mesh.half_length,vacuum_cells+1)
        zt=np.linspace(mesh.half_length,zout,vacuum_cells+1)
        self.rnodes=np.r_[mesh.rf,rs[1:]]; self.znodes=np.r_[zb,mesh.zf[1:],zt[1:]]
        nr=len(self.rnodes); nz=len(self.znodes); lo=vacuum_cells; hi=lo+mesh.nz
        # R wall, lower end, upper end; exclude axis and duplicated corners.
        self.bi=np.r_[np.full(mesh.nz+1,mesh.nr),np.arange(1,mesh.nr),np.arange(1,mesh.nr)]
        self.bj=np.r_[np.arange(lo,hi+1),np.full(mesh.nr-1,lo),np.full(mesh.nr-1,hi)]
        self.boundary_r=self.rnodes[self.bi]; self.boundary_z=self.znodes[self.bj]
        self.boundary_count=len(self.bi)
        ids=np.full((nr,nz),-1,int); ids[self.bi,self.bj]=np.arange(self.boundary_count)
        next_id=self.boundary_count
        for i in range(1,nr-1):
            for j in range(1,nz-1):
                if ids[i,j]>=0 or (i<mesh.nr and lo<j<hi):continue
                ids[i,j]=next_id;next_id+=1
        self.mode_start=next_id; self.size=next_id+modes; self.ids=ids
        self.patches=[]
        radial=(np.arange(mesh.nr+1),np.arange(mesh.nr,nr))
        axial=(np.arange(lo,hi+1),np.arange(0,lo+1),np.arange(hi,nz))
        for ri,ir in enumerate(radial):
            for zi,iz in enumerate(axial):
                if ri==0 and zi==0:continue
                side=(0,-1,1)[zi]
                r=self.rnodes[ir];z=self.znodes[iz];rr,zz=np.meshgrid(r,z,indexing='ij')
                patch_ids=ids[np.ix_(ir,iz)].ravel(); rows=[];cols=[];data=[]
                for k,node in enumerate(patch_ids):
                    if node>=0:rows.append(k);cols.append(node);data.append(1.)
                    elif rr.ravel()[k]!=0:
                        trace=multipoles(rr.ravel()[k],zz.ravel()[k],self.scale,modes)[0]
                        rows.extend([k]*modes);cols.extend(range(next_id,self.size));data.extend(trace)
                lift=csr_matrix((data,(rows,cols)),shape=(len(patch_ids),self.size))
                self.patches.append(dict(r=r,z=z,radial_side=ri,axial_side=side,lift=lift))
        matrix=np.zeros((self.size,self.size))
        for patch in self.patches:matrix+=self._energy(patch)
        matrix[next_id:,next_id:]+=2*np.pi/MU0*exterior_energy(self.outer,modes)
        skew=np.linalg.norm(matrix-matrix.T)/np.linalg.norm(matrix)
        if skew>1e-11:raise RuntimeError('Nonreciprocal finite vacuum energy assembly')
        self.full_hessian=.5*(matrix+matrix.T)
        nb=self.boundary_count; vv=self.full_hessian[nb:,nb:]; scale=1/np.sqrt(np.diag(vv))
        factor=cho_factor(scale[:,None]*vv*scale[None,:],lower=True)
        self.response=-scale[:,None]*cho_solve(factor,scale[:,None]*self.full_hessian[nb:,:nb])
        reduced=self.full_hessian[:nb,:nb]+self.full_hessian[:nb,nb:]@self.response
        self.hessian=.5*(reduced+reduced.T)
        if np.linalg.eigvalsh(self.hessian)[0]<=0:raise RuntimeError('Nonpositive exterior magnetic energy')
        self.reciprocity_relative=float(skew)

    def _trace_error(self, points, nodes, *, radial, side=1):
        """Multipole trace minus its nodal polynomial; derivatives share the map."""
        if radial:
            z=np.full_like(points,side*self.outer.half_length)
            value,dr,_=multipoles(points,z,self.scale,self.modes)
            nodal=multipoles(nodes,np.full_like(nodes,side*self.outer.half_length),self.scale,self.modes)[0]
            n,d,_,_=bases(.5*nodes**2,.5*points**2,self.degree)
            derivative=dr-points[:,None]*(d@nodal)
            ell=np.arange(1,self.modes+1);distance=abs(z)
            axis=(self.scale/distance[:,None])**ell/distance[:,None]**2*np.sqrt(ell*(ell+1))*np.sign(z[:,None])**(ell-1)
            bz=np.divide(dr,points[:,None],out=axis.copy(),where=points[:,None]>0)-d@nodal
        else:
            r=np.full_like(points,self.outer.radius)
            value,_,dz=multipoles(r,points,self.scale,self.modes)
            nodal=multipoles(np.full_like(nodes,self.outer.radius),nodes,self.scale,self.modes)[0]
            n,d,_,_=bases(nodes,points,self.degree)
            derivative=dz-d@nodal;bz=None
        return value-n@nodal,derivative,bz

    def _correction(self, patch, r, z, *, grid):
        shape=(len(r),len(z),self.modes) if grid else (len(r),self.modes)
        value=np.zeros(shape);dr=np.zeros(shape);dz=np.zeros(shape);bz_axis=np.zeros(shape)
        if patch['radial_side']:
            error,derivative,_=self._trace_error(z,patch['z'],radial=False)
            width=self.outer.radius-self.mesh.radius;blend=(r-self.mesh.radius)/width
            if grid:
                value+=blend[:,None,None]*error[None,:,:];dr+=error[None,:,:]/width
                dz+=blend[:,None,None]*derivative[None,:,:]
            else:
                value+=blend[:,None]*error;dr+=error/width;dz+=blend[:,None]*derivative
        side=patch['axial_side']
        if side:
            error,derivative,bz=self._trace_error(r,patch['r'],radial=True,side=side)
            width=self.outer.half_length-self.mesh.half_length;blend=(side*z-self.mesh.half_length)/width
            if grid:
                value+=error[:,None,:]*blend[None,:,None];dr+=derivative[:,None,:]*blend[None,:,None]
                dz+=error[:,None,:]*(side/width);bz_axis+=bz[:,None,:]*blend[None,:,None]
            else:
                value+=error*blend[:,None];dr+=derivative*blend[:,None];dz+=error*(side/width);bz_axis+=bz*blend[:,None]
        return value,dr,dz,bz_axis

    def _energy(self, patch):
        q,w=np.polynomial.legendre.leggauss(self.order)
        def quadrature(nodes):
            width=np.diff(nodes);return (nodes[:-1,None]+.5*width[:,None]*(q+1)).ravel(),(.5*width[:,None]*w).ravel()
        r,wr=quadrature(patch['r']);z,wz=quadrature(patch['z'])
        nr,dx,_,_=bases(.5*patch['r']**2,.5*r*r,self.degree)
        nz,dz,_,_=bases(patch['z'],z,self.degree)
        radial=kron(diags(r)@dx,nz,format='csr');axial=kron(nr,dz,format='csr')
        weight=(2*np.pi/MU0*(wr/r)[:,None]*wz[None,:]).ravel()
        _,cr,cz,_=self._correction(patch,r,z,grid=True);cr=cr.reshape(-1,self.modes);cz=cz.reshape(-1,self.modes)
        local=radial.T@radial.multiply(weight[:,None])+axial.T@axial.multiply(weight[:,None])
        lift=patch['lift'];result=(lift.T@local@lift).toarray()
        cross=lift.T@(radial.T@(weight[:,None]*cr)+axial.T@(weight[:,None]*cz))
        start=self.mode_start;result[:,start:]+=cross;result[start:,:]+=cross.T
        result[start:,start:]+=cr.T@(weight[:,None]*cr)+cz.T@(weight[:,None]*cz)
        return result

    def coefficients(self, boundary):
        boundary=np.asarray(boundary,float)
        if boundary.shape!=(self.boundary_count,) or not np.isfinite(boundary).all():
            raise ValueError('Finite complete cylinder boundary flux required')
        return np.r_[boundary,self.response@boundary]

    def samples(self, boundary, r, z):
        """Flux and Br,Bz at arbitrary exterior points, including the axis."""
        r,z=np.broadcast_arrays(np.asarray(r,float),np.asarray(z,float));shape=r.shape;r=r.ravel();z=z.ravel()
        if not np.isfinite(r).all() or not np.isfinite(z).all() or np.any(r<0):raise ValueError('Finite r>=0 exterior points required')
        if np.any((r<self.mesh.radius)&(abs(z)<self.mesh.half_length)):
            raise ValueError('Vacuum samples must be outside or on the finite cylinder')
        coefficients=self.coefficients(boundary);mode=coefficients[self.mode_start:]
        flux=np.zeros_like(r);br=np.zeros_like(r);bz=np.zeros_like(r)
        outer=(r>=self.outer.radius)|(abs(z)>=self.outer.half_length)
        if np.any(outer):
            rr,zz=r[outer],z[outer];p,pr,pz=multipoles(rr,zz,self.scale,self.modes)
            flux[outer]=p@mode;br[outer]=np.divide(-pz@mode,rr,out=np.zeros_like(rr),where=rr>0)
            ell=np.arange(1,self.modes+1);distance=np.hypot(rr,zz)
            axis=((self.scale/distance[:,None])**ell/distance[:,None]**2*np.sqrt(ell*(ell+1))*np.sign(zz[:,None])**(ell-1))@mode
            bz[outer]=np.divide(pr@mode,rr,out=axis,where=rr>0)
        rside=(r>=self.mesh.radius).astype(int);zside=np.where(z>=self.mesh.half_length,1,np.where(z<=-self.mesh.half_length,-1,0))
        for patch in self.patches:
            choose=(~outer)&(rside==patch['radial_side'])&(zside==patch['axial_side'])
            if not np.any(choose):continue
            rr,zz=r[choose],z[choose]
            nr,dx,_,_=bases(.5*patch['r']**2,.5*rr*rr,self.degree)
            nz,dz,_,_=bases(patch['z'],zz,self.degree)
            values=np.asarray(patch['lift']@coefficients).reshape(len(patch['r']),len(patch['z']))
            cp,cr,cz,cbz=self._correction(patch,rr,zz,grid=False)
            flux[choose]=np.sum((nr@values)*nz.toarray(),axis=1)+cp@mode
            pr=rr*np.sum((dx@values)*nz.toarray(),axis=1)+cr@mode
            pz=np.sum((nr@values)*dz.toarray(),axis=1)+cz@mode
            br[choose]=np.divide(-pz,rr,out=np.zeros_like(rr),where=rr>0)
            axis=np.sum((dx@values)*nz.toarray(),axis=1)+cbz@mode
            bz[choose]=np.divide(pr,rr,out=axis,where=rr>0)
        return flux.reshape(shape),br.reshape(shape),bz.reshape(shape)
