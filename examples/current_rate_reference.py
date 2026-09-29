"""Independent Fourier/Chebyshev inversion of the continuum toroidal inertia."""
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME,mu_0 as MU0,atomic_mass
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import LinearOperator,gmres
from .current_inertia_reference import fields,g_rate


class CurrentRateReference:
    """Regular analytic functions of s=r²/R², periodic z, physical g(R)=0."""
    def __init__(self,reference,radial_order=64,axial_points=128,ion_mass=39.948*atomic_mass):
        self.ref=reference;self.order=radial_order;self.nz=axial_points
        j=np.arange(radial_order+1);x=np.cos(np.pi*j/radial_order);self.s=(1-x)/2
        c=(-1.)**j;c[[0,-1]]*=2
        difference=x[:,None]-x[None,:]
        derivative=(c[:,None]/c[None,:])/(difference+np.eye(len(x)))
        derivative-=np.diag(derivative.sum(axis=1));self.derivative=-2*derivative
        self.second=self.derivative@self.derivative
        self.z=-reference.length/2+np.arange(axial_points)*reference.length/axial_points
        self.k=2*np.pi*np.fft.fftfreq(axial_points,d=reference.length/axial_points)
        self.mu=ME*ion_mass/(ME+ion_mass);self.delta=(ion_mass-ME)/(ion_mass+ME)
        self.d2=self.mu/(MU0*reference.density*E*E)
        s,z=self.s[:,None],self.z[None];k0=2*np.pi/reference.length
        self.relative_density=np.exp(reference.density_contrast*(1-s)**2*np.cos(k0*z))
        self.ascalar=-2*reference.density_contrast*(1-s)*np.cos(k0*z)
        self.az=-reference.density_contrast*k0*(1-s)**2*np.sin(k0*z)
        radial=-4*self.d2/reference.radius**2*(s*self.second+2*self.derivative)
        self.factors={}
        for k in self.k:
            if abs(k) not in self.factors:
                matrix=radial+np.eye(len(x))*(1+self.d2*k*k)
                matrix[-1]=0;matrix[-1,-1]=1
                self.factors[abs(k)]=lu_factor(matrix)

    def apply(self,value):
        s=self.s[:,None];ds=self.derivative@value;dss=self.second@value
        spectral=np.fft.fft(value,axis=1)
        dz=np.fft.ifft(1j*self.k[None]*spectral,axis=1).real
        dzz=np.fft.ifft(-self.k[None]**2*spectral,axis=1).real
        elliptic=-4/self.ref.radius**2*(s*dss+2*ds-self.ascalar*(s*ds+value))-dzz+self.az*dz
        result=value+self.d2/self.relative_density*elliptic;result[-1]=value[-1]
        return result

    def solve(self,source,tolerance=1e-11):
        shape=(self.order+1,self.nz);size=np.prod(shape);rhs=np.asarray(source,float).copy()
        if rhs.shape!=shape or not np.isfinite(rhs).all():raise ValueError('Finite collocation source required')
        rhs[-1]=0.
        def action(v):return self.apply(v.reshape(shape)).ravel()
        def precondition(v):
            transformed=np.fft.fft(v.reshape(shape),axis=1)
            solved=np.column_stack([lu_solve(self.factors[abs(k)],transformed[:,j]) for j,k in enumerate(self.k)])
            return np.fft.ifft(solved,axis=1).real.ravel()
        iterations=[]
        solved,info=gmres(LinearOperator((size,size),matvec=action),rhs.ravel(),
            M=LinearOperator((size,size),matvec=precondition),rtol=tolerance,atol=0.,restart=60,maxiter=20,
            callback=iterations.append,callback_type='pr_norm')
        residual=float(np.linalg.norm(action(solved)-rhs.ravel())/max(np.linalg.norm(rhs),np.finfo(float).tiny))
        if info or residual>2*tolerance:raise RuntimeError(f'Independent continuum inertia solve failed: {info}, {residual}')
        return solved.reshape(shape),dict(iterations=len(iterations),relative_residual=residual,tolerance=tolerance)

    def original_source(self,axis_probe_fraction=1e-8):
        r=self.ref.radius*np.sqrt(self.s[:,None]);r[0]=self.ref.radius*axis_probe_fraction;z=self.z[None]
        ref=self.ref;mu=self.mu;delta=self.delta
        def hall(rr,zz):
            w=fields(ref,rr,zz,mu,delta)['drift'];b=rr*ref.fields(rr,zz)['g']
            return delta*np.stack((-w[...,2]*b,np.zeros_like(b),w[...,0]*b),axis=-1)
        epsilon=1e-20
        hall_rate=(hall(r+1j*epsilon,z)[...,2].imag-hall(r,z+1j*epsilon)[...,0].imag)/(epsilon*r)
        beta=.5*(1+delta)
        return hall_rate+beta*ref.fields(r,z)['g_rate']+g_rate(ref,r,z,mu,delta)

    def cell_average(self,mesh,values,order=12):
        """Independent radial quadrature and exact Fourier cell integrals."""
        nodes,weights=np.polynomial.legendre.leggauss(order);wi=(-1.)**np.arange(self.order+1);wi[[0,-1]]*=.5
        lo=(mesh.rf[:-1]/self.ref.radius)**2;hi=(mesh.rf[1:]/self.ref.radius)**2
        radial=np.zeros((mesh.nr,self.order+1))
        for a,wa in zip(nodes,weights):
            point=.5*(lo+hi+a*(hi-lo));difference=point[:,None]-self.s[None]
            exact=difference==0;safe=np.where(exact,1.,difference);w=wi[None]/safe
            has=exact.any(axis=1);w[has]=exact[has]
            radial+=.5*wa*w/w.sum(axis=1)[:,None]
        axial=np.exp(1j*(mesh.z[:,None]+self.ref.length/2)*self.k[None])*np.sinc(self.k[None]*mesh.dz/(2*np.pi))
        return (radial@np.fft.fft(values,axis=1)@axial.T).real/self.nz
