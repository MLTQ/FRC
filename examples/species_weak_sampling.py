"""Independent tensor weak/face actions for the species entropy budget."""
import numpy as np
from numpy.polynomial import Polynomial
from .current_transport_sampling import TransportSampler


def fourier(fluid, points):
    n=fluid.nz; frequency=np.fft.fftfreq(n)*n
    phase=(np.asarray(points)+fluid.length/2)/fluid.length
    values=np.exp(2j*np.pi*phase[:,None]*frequency)
    derivative=values*(2j*np.pi*frequency/fluid.length)
    if n%2==0:
        values[:,n//2]=np.cos(np.pi*n*phase)
        derivative[:,n//2]=-np.pi*n/fluid.length*np.sin(np.pi*n*phase)
    coefficients=np.fft.fft(np.eye(n),axis=1)/n
    return (coefficients@values.T).real.T,(coefficients@derivative.T).real.T


def radial_faces(fluid, mesh):
    radius=np.column_stack((mesh.rf[:-1],mesh.rf[1:])).ravel()
    owner=np.repeat(np.arange(mesh.nr)//(mesh.nr//fluid.elements),2)
    result=np.zeros((len(radius),fluid.elements*(fluid.degree+1)))
    for row,(r,e) in enumerate(zip(radius,owner)):
        nodes=fluid.x[e,:,0];scale=nodes[-1]-nodes[0]
        roots=(nodes-nodes[0])/scale;point=(r*r/2-nodes[0])/scale
        for i in range(len(nodes)):
            others=np.delete(roots,i)
            result[row,e*(fluid.degree+1)+i]=Polynomial.fromroots(others)(point)/np.prod(roots[i]-others)
    return radius,result


class WeakSampler:
    def __init__(self,fluid,mesh,order):
        self.sample=TransportSampler(fluid,mesh,order)
        self.fluid=fluid;self.mesh=mesh;self.order=order
        s=self.sample;self.R,self.DR=s.radial[False];self.Z,self.DZ=fourier(fluid,s.z)
        self.NR=s.nx;self.NDR=s.R*s.dx;self.NZ=s.nz;self.NDZ=s.dz
        _,w=np.polynomial.legendre.leggauss(order)
        self.wr=(np.pi*mesh.dr*s.r.reshape(mesh.nr,order)*w).ravel()
        self.wz=np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()
        self.projection=(np.linalg.solve(self.NR.T@(self.wr[:,None]*self.NR),
                                        self.NR.T@(self.wr[:,None]*self.R)),
                         np.linalg.solve(self.NZ.T@(self.wz[:,None]*self.NZ),
                                         self.NZ.T@(self.wz[:,None]*self.Z)))
        self.rf,self.RF=radial_faces(fluid,mesh)
        self.zf=np.column_stack((mesh.zf[:-1],mesh.zf[1:])).ravel()
        self.ZF,_=fourier(fluid,self.zf)

    def shape(self,value):
        return np.asarray(value).reshape(self.R.shape[1],self.fluid.nz)

    def scalar(self,value):
        return self.R@self.shape(value)@self.Z.T

    def adjoint(self,value):
        return self.R.T@(self.sample.weight*value)@self.Z

    def gradient_adjoint(self,radial,axial,projection=None):
        if projection is None:
            return self.DR.T@(self.sample.weight*radial)@self.Z+self.R.T@(self.sample.weight*axial)@self.DZ
        pr,pz=projection
        potential=self.NDR.T@(self.sample.weight*radial)@self.NZ+self.NR.T@(self.sample.weight*axial)@self.NDZ
        return pr.T@potential@pz

    def faces(self,entropy,radial_current,axial_current):
        """Every one-sided outward covector; no assumptions about cancellations."""
        entropy=self.shape(entropy)
        jr=radial_current.reshape(len(self.rf),len(self.sample.z))
        jz=axial_current.reshape(len(self.sample.r),len(self.zf))
        radial=(self.RF@entropy@self.Z.T)*jr*(2*np.pi*self.rf[:,None]*self.wz)
        radial*=np.tile([-1.,1.],self.mesh.nr)[:,None]
        axial=(self.R@entropy@self.ZF.T)*jz*self.wr[:,None]
        axial*=np.tile([-1.,1.],self.mesh.nz)[None,:]
        out={}
        for name,ids in [('radial_axis',[0]),('radial_wall',[len(self.rf)-1]),
                         ('radial_interior',np.arange(1,len(self.rf)-1))]:
            out[name]=self.RF[ids].T@radial[ids]@self.Z
        for name,ids in [('axial_seam',[0,len(self.zf)-1]),
                         ('axial_interior',np.arange(1,len(self.zf)-1))]:
            out[name]=self.R.T@axial[:,ids]@self.ZF[ids]
        return out

    def decompose(self,entropy,current,divergence,radial_current,axial_current,projection):
        s=self.shape(entropy);sh=self.scalar(s)
        sr=self.DR@s@self.Z.T;sz=self.R@s@self.DZ.T
        pr,pz=projection;potential=pr@s@pz.T
        phr=self.NDR@potential@self.NZ.T;phz=self.NR@potential@self.NDZ.T
        jr,jz=current[...,0],current[...,2]
        reverse=self.gradient_adjoint(jr*sh,jz*sh,projection)
        direct=.5*(self.adjoint(jr*phr+jz*phz)-reverse)
        terms={'advection':self.adjoint(jr*sr+jz*sz),
               'volume_charge':.5*self.adjoint(sh*divergence),
               'transfer':.5*(self.adjoint(jr*(phr-sr)+jz*(phz-sz))-reverse+
                              self.gradient_adjoint(jr*sh,jz*sh))}
        terms.update({name:-.5*value for name,value in self.faces(s,radial_current,axial_current).items()})
        terms['sum']=sum(terms.values());terms['direct']=direct
        return terms
