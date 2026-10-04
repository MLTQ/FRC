"""Independent root-polynomial/Fourier sampler for the transport audit."""
import numpy as np
from numpy.polynomial import Polynomial


def bases(edges, points, degree, periodic=False):
    count=len(edges)-1;size=count if periodic else count+1
    nodal=[np.zeros((len(points),size)) for _ in range(3)]
    edge=[np.zeros((len(points),count)) for _ in range(2)]
    for row,point in enumerate(points):
        if periodic:point=edges[0]+(point-edges[0])%(edges[-1]-edges[0])
        cell=min(np.searchsorted(edges,point,side='right')-1,count-1)
        start=cell-(degree-1)//2
        if not periodic:start=min(max(start,0),count-degree)
        ids=np.arange(start,start+degree+1)
        coords=edges[0]+ids*(edges[1]-edges[0]) if periodic else edges[ids]
        scale=coords[-1]-coords[0];roots=(coords-point)/scale;derivatives=[]
        for i,index in enumerate(ids):
            others=np.delete(roots,i);poly=Polynomial.fromroots(others)/np.prod(roots[i]-others)
            values=[poly.deriv(k)(0.)/scale**k for k in range(3)]
            for k in range(3):nodal[k][row,index%size]+=values[k]
            derivatives.append(values)
        for i,index in enumerate(ids[:-1]):
            for k in range(2):edge[k][row,index%count]=-sum(v[k+1] for v in derivatives[:i+1])*(coords[i+1]-coords[i])
    return *nodal,*edge


class TransportSampler:
    def __init__(self,fluid,mesh,order):
        self.d=fluid;self.mesh=mesh;q,w=np.polynomial.legendre.leggauss(order)
        self.r=(mesh.r[:,None]+.5*mesh.dr*q).ravel();self.z=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
        self.R=self.r[:,None];self.Z=self.z[None,:]
        self.weight=(np.pi*mesh.dr*self.r.reshape(mesh.nr,-1)*w).ravel()[:,None]*np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()[None,:]
        self.nx,self.dx,_,self.ex,self.dex=bases(.5*mesh.rf**2,.5*self.r**2,5)
        self.nz,self.dz,_,self.ez,self.dez=bases(mesh.zf,self.z,5,True)
        d=fluid;phase=(self.z+.5*d.length)/d.length;freq=np.fft.fftfreq(d.nz)*d.nz
        self.fourier=np.exp(2j*np.pi*phase[:,None]*freq[None,:]);self.fourier_d=self.fourier*(2j*np.pi*freq/d.length)
        if d.nz%2==0:
            self.fourier[:,d.nz//2]=np.cos(np.pi*d.nz*phase)
            self.fourier_d[:,d.nz//2]=-np.pi*d.nz/d.length*np.sin(np.pi*d.nz*phase)
        self.radial={};element=np.minimum((self.r*d.elements/d.radius).astype(int),d.elements-1)
        for transverse in (False,True):
            value=np.zeros((len(self.r),d.elements*(d.degree+1)));derivative=value.copy()
            for e in range(d.elements):
                rows=np.flatnonzero(element==e);begin=1 if transverse and e==0 else 0
                nodes=d.x[e,begin:,0];scale=nodes[-1]-nodes[0];origin=nodes[0];roots=(nodes-origin)/scale;points=(.5*self.r[rows]**2-origin)/scale
                for i in range(len(nodes)):
                    others=np.delete(roots,i);poly=Polynomial.fromroots(others)/np.prod(roots[i]-others);col=e*(d.degree+1)+begin+i
                    value[rows,col]=poly(points);derivative[rows,col]=poly.deriv()(points)*self.r[rows]/scale
            self.radial[transverse]=(value,derivative)

    def scalar(self,value,transverse=False,axis=None):
        matrix=self.radial[transverse][axis==0];radial=matrix@np.asarray(value).reshape(-1,self.d.nz)
        result=(np.fft.fft(radial,axis=1)/self.d.nz)@(self.fourier_d if axis==2 else self.fourier).T
        if np.max(abs(result.imag))>1e-10*max(np.max(abs(result.real)),np.finfo(float).tiny):raise ValueError('Nonreal independent scalar')
        return result.real

    def electric(self,value,gradient=False):
        m=self.mesh;r=self.R;count=m.nr*m.nz;theta=(m.nr-1)*m.nz
        jr=value[:count].reshape(m.nr,m.nz)
        jt=np.pad(value[count:count+theta].reshape(m.nr-1,m.nz),((1,1),(0,0)))*m.rf[:,None]
        jz=value[count+theta:].reshape(m.nr+1,m.nz)
        out=np.stack((r*(self.ex@jr@self.nz.T),(self.nx@jt@self.nz.T)/r,self.nx@jz@self.ez.T),axis=-1)
        if not gradient:return out
        grad=np.zeros(out.shape+(3,))
        grad[...,0]=np.stack(((self.ex+r*r*self.dex)@jr@self.nz.T,(self.dx-self.nx/(r*r))@jt@self.nz.T,r*(self.dx@jz@self.ez.T)),axis=-1)
        grad[...,2]=np.stack((r*(self.ex@jr@self.dz.T),(self.nx@jt@self.dz.T)/r,self.nx@jz@self.dez.T),axis=-1)
        grad[...,0,1]=-out[...,1]/r;grad[...,1,1]=out[...,0]/r
        return out,grad

    def magnetic(self,value):
        m=self.mesh;count=(m.nr+1)*m.nz;p=value[:count].reshape(m.nr+1,m.nz);g=value[count:].reshape(m.nr,m.nz)
        return np.stack((-(self.nx@p@self.dz.T)/self.R,self.R*(self.ex@g@self.ez.T),self.dx@p@self.nz.T),axis=-1)

    def fluid(self,state,rate=None):
        n=self.scalar(state[...,0]);r=self.R
        if np.any(n<=0):raise ValueError('Nonpositive independently sampled density')
        v=np.stack((r*self.scalar(state[...,1]/state[...,0],True),r*self.scalar(state[...,2]/state[...,0],True),self.scalar(state[...,3])/n),axis=-1)
        grad=np.zeros(v.shape+(3,))
        for axis in (0,2):
            for k in (0,1):grad[...,k,axis]=r*self.scalar(state[...,k+1]/state[...,0],True,axis)+(v[...,k]/r if axis==0 else 0)
            grad[...,2,axis]=(self.scalar(state[...,3],axis=axis)-v[...,2]*self.scalar(state[...,0],axis=axis))/n
        grad[...,0,1]=-v[...,1]/r;grad[...,1,1]=v[...,0]/r
        if rate is None:return n,v,grad
        nt=self.scalar(rate[...,0]);vt=np.empty_like(v)
        for k in (0,1):vt[...,k]=r*self.scalar((rate[...,k+1]-state[...,k+1]*rate[...,0]/state[...,0])/state[...,0],True)
        vt[...,2]=(self.scalar(rate[...,3])-v[...,2]*nt)/n
        return n,v,grad,nt,vt
