"""Shape-dependent magnetic-only exterior for the free material cylinder."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.linalg import cho_factor,cho_solve
from scipy.sparse import csr_matrix,diags,hstack,kron,vstack
from .compatible_polynomials import bases
from .finite_cylinder_vacuum import FiniteCylinderVacuum
from .open_boundary import exterior_energy


class MaterialMagneticVacuum(FiniteCylinderVacuum):
    """Condense vacuum flux on an explicitly deformed, massless exterior mesh.

    Boundary coefficients are psi, using the parent's polynomial trace. They
    must not be spliced onto a different r0²*a interpolation at the plasma ends.
    """
    def __init__(self,mesh,**options):
        super().__init__(mesh,**options)
        rr,zz=mesh.corners;self.reference=np.column_stack((rr.ravel(),zz.ravel()))
        self.boundary_nodes=np.r_[np.arange(mesh.nr*(mesh.nz+1),(mesh.nr+1)*(mesh.nz+1)),
            np.arange(1,mesh.nr)*(mesh.nz+1),np.arange(1,mesh.nr)*(mesh.nz+1)+mesh.nz]
        self.outer_hessian=np.zeros_like(self.full_hessian)
        self.outer_hessian[self.mode_start:,self.mode_start:]=2*np.pi/MU0*exterior_energy(self.outer,self.modes)
        q,w=np.polynomial.legendre.leggauss(self.order)
        def quadrature(nodes):
            width=np.diff(nodes)
            return (nodes[:-1,None]+.5*width[:,None]*(q+1)).ravel(),(.5*width[:,None]*w).ravel()
        points=[];weights=[];values=[];radials=[];axials=[];flux_r=[];flux_z=[]
        for patch in self.patches:
            rn,zn=patch['r'],patch['z'];r,wr=quadrature(rn);z,wz=quadrature(zn)
            rgrid,zgrid=np.meshgrid(r,z,indexing='ij');points.append(np.column_stack((rgrid.ravel(),zgrid.ravel())))
            weights.append((wr[:,None]*wz[None,:]).ravel())
            nr,dr,_,_=bases(rn,r,1);nz,dz,_,_=bases(zn,z,1)
            lift=self._geometry_lift(rn,zn)
            values.append(kron(nr,nz,format='csr')@lift)
            radials.append(kron(dr,nz,format='csr')@lift);axials.append(kron(nr,dz,format='csr')@lift)
            nr,dx,_,_=bases(.5*rn**2,.5*r*r,self.degree);nz,dz,_,_=bases(zn,z,self.degree)
            _,cr,cz,_=self._correction(patch,r,z,grid=True)
            zeros=csr_matrix((len(r)*len(z),self.mode_start))
            flux_r.append(kron(diags(r)@dx,nz,format='csr')@patch['lift']+hstack((zeros,csr_matrix(cr.reshape(-1,self.modes)))))
            flux_z.append(kron(nr,dz,format='csr')@patch['lift']+hstack((zeros,csr_matrix(cz.reshape(-1,self.modes)))))
        self.points0=np.concatenate(points);self.weights0=np.concatenate(weights)
        self.motion=vstack(values,format='csr');self.motion_r=vstack(radials,format='csr');self.motion_z=vstack(axials,format='csr')
        self.flux_r=vstack(flux_r,format='csr');self.flux_z=vstack(flux_z,format='csr')

    def _geometry_lift(self,r,z):
        """Q1 boundary motion extends to zero at the fixed outer matching box."""
        mesh=self.mesh;rows=[];cols=[];values=[]
        for i,ri in enumerate(r):
            outside_r=ri>=mesh.radius
            ir=mesh.nr if outside_r else int(np.argmin(abs(mesh.rf-ri)))
            ar=(self.outer.radius-ri)/(self.outer.radius-mesh.radius) if outside_r else 1.
            for j,zj in enumerate(z):
                outside_z=abs(zj)>=mesh.half_length
                iz=(0 if zj<0 else mesh.nz) if outside_z else int(np.argmin(abs(mesh.zf-zj)))
                az=(self.outer.half_length-abs(zj))/(self.outer.half_length-mesh.half_length) if outside_z else 1.
                if ar*az:
                    rows.append(i*len(z)+j);cols.append(ir*(mesh.nz+1)+iz);values.append(ar*az)
        return csr_matrix((values,(rows,cols)),shape=(len(r)*len(z),len(self.reference)))

    def geometry(self,coordinates,velocity=None):
        coordinates=np.asarray(coordinates,float)
        if coordinates.shape!=self.reference.shape or not np.isfinite(coordinates).all():raise ValueError('Complete finite material coordinates required')
        if np.any(coordinates[self.reference[:,0]==0,0]!=0):raise ValueError('Axis radial displacement must remain zero')
        displacement=coordinates-self.reference;points=self.points0+self.motion@displacement;r=points[:,0]
        F=np.empty((len(r),2,2));F[:,:,0]=self.motion_r@displacement;F[:,:,1]=self.motion_z@displacement
        F+=np.eye(2);det=np.linalg.det(F)
        if np.any(det<=0) or np.any(r<=0):raise ValueError('Moving vacuum mesh folded or crossed the axis; no mesh/state repair')
        inv=np.linalg.inv(F)
        gr=diags(inv[:,0,0])@self.flux_r+diags(inv[:,1,0])@self.flux_z
        gz=diags(inv[:,0,1])@self.flux_r+diags(inv[:,1,1])@self.flux_z
        nr=diags(inv[:,0,0])@self.motion_r+diags(inv[:,1,0])@self.motion_z
        nz=diags(inv[:,0,1])@self.motion_r+diags(inv[:,1,1])@self.motion_z
        weight=2*np.pi/MU0*self.weights0*det/r
        full=(gr.T@gr.multiply(weight[:,None])+gz.T@gz.multiply(weight[:,None])).toarray()+self.outer_hessian
        nb=self.boundary_count;vv=full[nb:,nb:];scale=1/np.sqrt(np.diag(vv))
        factor=cho_factor(scale[:,None]*vv*scale[None,:],lower=True)
        response=-scale[:,None]*cho_solve(factor,scale[:,None]*full[nb:,:nb])
        extension=np.vstack((np.eye(nb),response));h=full[:nb,:nb]+full[:nb,nb:]@response
        result=dict(points=points,r=r,determinant=det,weight=weight,gradient_r=gr,gradient_z=gz,
            motion_r=nr,motion_z=nz,full_hessian=full,extension=extension,hessian=.5*(h+h.T))
        if velocity is not None:
            velocity=np.asarray(velocity,float)
            if velocity.shape!=coordinates.shape or not np.isfinite(velocity).all():raise ValueError('Complete finite boundary velocity required')
            grad=np.stack((nr@velocity,nz@velocity),axis=-1);vr=(self.motion@velocity)[:,0]
            diagonal=grad[:,0,0]+grad[:,1,1]-vr/r
            metric=diagonal[:,None,None]*np.eye(2)-grad-grad.transpose(0,2,1)
            operators=(gr,gz)
            full_rate=sum(operators[i].T@operators[j].multiply((weight*metric[:,i,j])[:,None]) for i in range(2) for j in range(2))
            rate=extension.T@(full_rate@extension)
            result.update(hessian_rate=.5*(rate+rate.T),velocity=velocity)
        return result

    def evaluate(self,coordinates,boundary,velocity=None):
        return self.field_response(self.geometry(coordinates,velocity),boundary)

    def field_response(self,g,boundary):
        """Energy and shape force from one already assembled actual geometry."""
        boundary=np.asarray(boundary,float)
        if boundary.shape!=(self.boundary_count,) or not np.isfinite(boundary).all():raise ValueError('Finite complete physical boundary flux required')
        g=dict(g);c=g['extension']@boundary
        pr=g['gradient_r']@c;pz=g['gradient_z']@c;weight=g['weight'];r=g['r'];square=pr*pr+pz*pz
        nr,nz=g['motion_r'],g['motion_z']
        radial=nr.T@(.5*weight*(pz*pz-pr*pr))-nz.T@(weight*pr*pz)-self.motion.T@(.5*weight*square/r)
        axial=nz.T@(.5*weight*(pr*pr-pz*pz))-nr.T@(weight*pr*pz)
        derivative=np.column_stack((radial,axial));derivative[self.reference[:,0]==0,0]=0.
        g.update(coefficients=c,magnetic_field=np.column_stack((-pz/r,pr/r)),
            energy=float(.5*boundary@g['hessian']@boundary),boundary_load=g['hessian']@boundary,
            coordinate_energy_gradient=derivative,force=-derivative,
            volume=2*np.pi*r*self.weights0*g['determinant'])
        if 'velocity' in g:g['geometry_power']=float(np.sum(derivative*g['velocity']))
        return g

    def samples(self,*args,**kwargs):
        raise ValueError('Use evaluate on the actual moving geometry; fixed-cylinder samples are inapplicable')
