"""One physical radial momentum field for continuity, velocity and kinetic work."""
import numpy as np
from scipy.sparse import kron
from scipy.sparse.linalg import LinearOperator,cg,splu
from .compatible_geometry import gram
from .compatible_velocity import radial_moments,axial_moments
from .particle_flux_basis import primitive_basis,difference
from .momentum_derivatives import radial_derivative,axial_derivative


class FaceRadialKinematics:
    """Interior radial faces store r n vr; cell momenta store (n vz,n r vtheta)."""
    def __init__(self,geometry,total_mass):
        if not np.isfinite(total_mass) or total_mass<=0:raise ValueError('Positive total species mass required')
        self.geometry=geometry;self.mesh=geometry.mesh;self.mass=total_mass;g=geometry;m=self.mesh
        _,hx=primitive_basis(.5*m.rf**2,.5*g.radial_points**2,g.degree)
        self.radial=hx[:,1:-1].multiply(1/g.radial_points[:,None]).tocsr()
        self.angular=radial_moments(m,g.radial_points,g.degree,power=2)
        self.axial=axial_moments(m,g.axial_points,g.degree);self.dx=difference(.5*m.rf**2)
        self.factor=splu(kron(gram(self.radial,g.wr),gram(self.axial,g.wz),format='csc'))
        self.dr={p:radial_derivative(m,g.radial_points,g.degree,p) for p in (0,2)}
        self.dz=axial_derivative(m,g.axial_points,g.degree)

    def validate(self,density,radial,momentum):
        m=self.mesh;q=np.asarray(radial,float);p=np.asarray(momentum,float)
        if density.geometry is not self.geometry:raise ValueError('One common density/kinetic geometry required')
        if q.shape!=(m.nr-1,m.nz) or p.shape!=(m.nr,m.nz,2) or not np.isfinite(q).all() or not np.isfinite(p).all():
            raise ValueError('Finite interior radial faces and two conservative cell momenta required')
        return q,p

    def full_faces(self,radial):return np.pad(radial,((1,1),(0,0)))

    def sampled_momentum(self,density,radial,momentum):
        q,p=self.validate(density,radial,momentum);g=self.geometry
        return np.stack((g.tensor(self.radial,q,density.axial),
                         g.tensor(self.angular,p[...,1],density.axial)/g.r,
                         g.tensor(density.radial,p[...,0],density.axial)),axis=-1)

    def fields(self,density,radial,momentum):
        p=self.sampled_momentum(density,radial,momentum);v=p/density.samples[...,None]
        energy=.5*self.mass*np.sum(p*v,axis=-1)
        if not np.isfinite(v).all() or not np.isfinite(energy).all():raise ValueError('Kinetic state exceeds finite range')
        return dict(momentum=p,velocity=v,kinetic=energy)

    def radial_adjoint(self,density,values):
        g=self.geometry
        return np.asarray(self.radial.T@g._unlayout(g.weight*values)@density.axial)

    def cell_adjoint(self,density,values):
        g=self.geometry
        axial=np.asarray(density.radial.T@g._unlayout(g.weight*values[...,2])@density.axial)
        angular=np.asarray(self.angular.T@g._unlayout(g.weight*values[...,1]/g.r)@density.axial)
        return np.stack((axial,angular),axis=-1)

    def energy_gradients(self,density,radial,momentum):
        v=self.fields(density,radial,momentum)['velocity']
        return dict(radial=self.mass*self.radial_adjoint(density,v[...,0]),
                    momentum=self.mass*self.cell_adjoint(density,v),
                    density=density.adjoint(-.5*self.mass*np.sum(v*v,axis=-1)))

    def average(self,value):return np.sum(self.geometry.weight*value,axis=(2,3))/self.mesh.volumes

    def energy_rate(self,density,radial,momentum,nt,qt,mt):
        v=self.fields(density,radial,momentum)['velocity'];pt=self.sampled_momentum(density,qt,mt)
        changing_momentum=self.mass*np.sum(v*pt,axis=-1)
        changing_density=-.5*self.mass*np.sum(v*v,axis=-1)*density.sample(nt)
        return dict(momentum=self.average(changing_momentum),density=self.average(changing_density),
                    total=self.average(changing_momentum+changing_density))

    def gradients(self,density,radial,momentum):
        q,p=self.validate(density,radial,momentum);g=self.geometry;n=density.samples
        f=self.fields(density,q,p);v=f['velocity']
        nr=g.tensor(self.dr[0],density.cell_means,density.axial);nz=g.tensor(density.radial,density.cell_means,self.dz)
        mr=g.tensor(density.radial,self.dx@self.full_faces(q),density.axial)-f['momentum'][...,0]/g.r
        pr=np.stack((mr,g.tensor(self.dr[2],p[...,1],density.axial)/g.r-f['momentum'][...,1]/g.r,
                     g.tensor(self.dr[0],p[...,0],density.axial)),axis=-1)
        pz=np.stack((g.tensor(self.radial,q,self.dz),g.tensor(self.angular,p[...,1],self.dz)/g.r,
                     g.tensor(density.radial,p[...,0],self.dz)),axis=-1)
        vr=(pr-v*nr[...,None])/n[...,None];vz=(pz-v*nz[...,None])/n[...,None]
        return dict(radial=vr,axial=vz,divergence=vr[...,0]+v[...,0]/g.r+vz[...,2],
                    curl=np.stack((-vz[...,1],vz[...,0]-vr[...,2],vr[...,1]+v[...,1]/g.r),axis=-1))

    def solve_rate(self,density,physical_advection,face_work):
        """Project physical particle-momentum advection minus generalized work."""
        g=self.geometry;shape=(self.mesh.nr-1,self.mesh.nz);n=density.samples
        scale=float(np.sum(g.weight*n)/self.mesh.volumes.sum());inverse=scale/n
        def action(flat):
            value=g.tensor(self.radial,np.asarray(flat).reshape(shape),density.axial)
            return self.radial_adjoint(density,inverse*value).ravel()
        load=self.radial_adjoint(density,inverse*physical_advection)-scale/self.mass*face_work[1:-1]
        if not np.isfinite(load).all():raise ValueError('Radial physical work exceeds finite range')
        if not np.any(load):return np.zeros(shape),0.
        op=LinearOperator((load.size,load.size),matvec=action,dtype=float)
        pre=LinearOperator(op.shape,matvec=self.factor.solve,dtype=float)
        value,info=cg(op,load.ravel(),M=pre,rtol=1e-13,atol=0.,maxiter=500)
        residual=float(np.linalg.norm(action(value)-load.ravel())/np.linalg.norm(load))
        if info or not np.isfinite(value).all() or residual>2e-12:
            raise RuntimeError(f'Physical radial face-momentum solve failed: {info}, {residual}')
        return value.reshape(shape),residual
