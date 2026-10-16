"""Linearized Landau heat transport in an orthonormal velocity basis."""
from math import factorial
import numpy as np
from scipy.linalg import eig, eigh, solve
from scipy.special import eval_hermitenorm, eval_laguerre, eval_genlaguerre, roots_hermitenorm, roots_laguerre


class KineticHeat:
    """Uniform Maxwellian, stationary ions, quasineutral parallel perturbations."""

    def __init__(self, degree, *, charge=1., collisions=True):
        if not isinstance(degree, int) or degree < 3 or not np.isfinite(charge) or charge <= 0:
            raise ValueError('Polynomial degree >=3 and positive ion charge required')
        self.degree=degree;self.charge=charge
        self.indices=[(d-2*q,q) for d in range(degree+1) for q in range(d//2+1)]
        self.count=len(self.indices);self.lookup={pair:i for i,pair in enumerate(self.indices)}
        self.density=self.lookup[(0,0)];self.momentum=self.lookup[(1,0)]
        self.energy=self.lookup[(2,0)];self.stress=self.lookup[(0,1)]
        rotation=np.eye(self.count);i,j=self.energy,self.stress
        rotation[np.ix_([i,j],[i,j])]=[[1/np.sqrt(3),np.sqrt(2/3)],[-np.sqrt(2/3),1/np.sqrt(3)]]
        streaming=np.zeros((self.count,self.count))
        for i,(p,q) in enumerate(self.indices):
            if (p+1,q) in self.lookup:
                j=self.lookup[(p+1,q)];streaming[i,j]=streaming[j,i]=np.sqrt(p+1)
        self.streaming=rotation.T@streaming@rotation
        self.heat=np.sqrt(1.5)*self.streaming[self.energy].copy()
        self.heat[self.momentum]=0. # Heat is energy flux minus 5/2 times particle flux.
        self.active=np.array([i for i in range(self.count) if i not in (self.density,self.momentum)])
        self.self_collision,self.ion_collision=self._collisions() if collisions else (np.zeros_like(streaming),np.zeros_like(streaming))

    def _gradient(self, points):
        x,y,z=points.T;r2=(x*x+y*y)/2
        h=[eval_hermitenorm(p,z)/np.sqrt(float(factorial(p))) for p in range(self.degree+1)]
        lag=[eval_laguerre(q,r2) for q in range(self.degree//2+1)]
        deriv=[np.zeros_like(x)]+[-eval_genlaguerre(q-1,1,r2) for q in range(1,self.degree//2+1)]
        result=np.empty((len(x),3,self.count))
        for i,(p,q) in enumerate(self.indices):
            result[:,0,i]=x*deriv[q]*h[p];result[:,1,i]=y*deriv[q]*h[p]
            result[:,2,i]=np.sqrt(p)*h[p-1]*lag[q] if p else 0.
        a=result[:,:,self.energy].copy();b=result[:,:,self.stress].copy()
        result[:,:,self.energy]=a/np.sqrt(3)-np.sqrt(2/3)*b
        result[:,:,self.stress]=np.sqrt(2/3)*a+b/np.sqrt(3)
        return result

    @staticmethod
    def _gram(gradient,direction,weight):
        parallel=np.einsum('nci,nc->ni',gradient,direction)
        transverse=gradient-direction[:,:,None]*parallel[:,None,:]
        result=np.zeros((gradient.shape[-1],gradient.shape[-1]))
        for axis in range(3):
            g=transverse[:,axis,:]*np.sqrt(weight)[:,None]
            result+=g.T@g
        return result

    def _collisions(self):
        d=self.degree;mu,wm=np.polynomial.legendre.leggauss(d+2)
        radial,wr=roots_laguerre((d+2)//2)
        # Rotational symmetry removes g azimuth; V still integrates all 3 axes.
        h,wh=roots_hermitenorm(d);h=h/np.sqrt(2);wh=wh/np.sqrt(2*np.pi)
        vx,vy,vz=np.meshgrid(h,h,h,indexing='ij');centre=np.column_stack([vx.ravel(),vy.ravel(),vz.ravel()])
        weights=(wh[:,None,None]*wh[None,:,None]*wh[None,None,:]).ravel()
        ee=np.zeros((self.count,self.count));ei=np.zeros_like(ee)
        for x,radial_weight in zip(radial,wr):
            for cosine,angular_weight in zip(mu,wm/2):
                direction=np.array([np.sqrt(1-cosine*cosine),0.,cosine])
                delta=np.sqrt(x)*direction
                gradient=self._gradient(centre+delta)-self._gradient(centre-delta)
                # Analytic collision invariants of the weak integrand.
                gradient[:,:,self.energy]=0.;gradient[:,:,self.momentum]=0.
                directions=np.broadcast_to(direction,centre.shape)
                ee+=.25*radial_weight*angular_weight/np.sqrt(np.pi)*self._gram(gradient,directions,weights)
                point=np.sqrt(2*x)*direction[None,:];gradient=self._gradient(point)
                gradient[:,:,self.energy]=0.
                ei+=self.charge/2*np.sqrt(2/np.pi)*radial_weight*angular_weight*self._gram(gradient,direction[None,:],np.ones(1))
        return ee,ei

    def local_conductivity(self, *, lorentz_only=False):
        odd=np.array([i for i,(p,q) in enumerate(self.indices) if p%2 and i!=self.momentum])
        collision=self.ion_collision if lorentz_only else self.ion_collision+self.self_collision
        heat=self.heat[odd]
        return float(heat@solve(collision[np.ix_(odd,odd)],heat,assume_a='sym'))

    def evolve(self, wave_number, times):
        """Return every requested Fourier-state sample, without damping/repair."""
        t=np.asarray(times,float);k=float(wave_number)
        if t.ndim!=1 or np.any(t<0) or not np.isfinite(t).all() or not np.isfinite(k):
            raise ValueError('Finite nonnegative sample times and finite K required')
        ix=np.ix_(self.active,self.active);c=(self.self_collision+self.ion_collision)[ix]
        s=self.streaming[ix];a0=np.zeros(len(self.active));a0[np.where(self.active==self.energy)[0][0]]=np.sqrt(1.5)
        if not np.any(c):
            frequencies,vectors=eigh(s);values=-1j*k*frequencies;amplitudes=vectors.T@a0;conditioning=1.
        else:
            values,vectors=eig(-c-1j*k*s);amplitudes=solve(vectors,a0);conditioning=float(np.linalg.cond(vectors))
        active=((np.exp(t[:,None]*values)*amplitudes)@vectors.T)
        states=np.zeros((len(t),self.count),complex);states[:,self.active]=active
        temperature=np.sqrt(2/3)*states[:,self.energy];heat=states@self.heat
        dissipation=2*np.einsum('ti,ij,tj->t',active.conj(),c,active).real
        return dict(time=t,state=states,temperature=temperature,heat=heat,
                    quadratic_free_energy=np.sum(abs(states)**2,axis=1),collision_dissipation=dissipation,
                    eigenvector_condition=conditioning,maximum_eigenvalue_real=float(values.real.max()))
