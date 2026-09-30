"""Independent root/FFT physical reconstruction for saved coupling arrays."""
import numpy as np
from numpy.polynomial import Polynomial
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .compatible_field_norms import direct_basis


def root_basis(nodes,points):
    scale=nodes[-1]-nodes[0];origin=nodes[0];x=(nodes-origin)/scale;y=(points-origin)/scale
    out=np.empty((len(points),len(nodes)))
    for i in range(len(nodes)):
        roots=np.delete(x,i);out[:,i]=Polynomial.fromroots(roots)(y)/np.prod(x[i]-roots)
    return out


def measure(d,mesh,ref,arrays,degree=5,order=16):
    q,w=np.polynomial.legendre.leggauss(order);r=(mesh.r[:,None]+.5*mesh.dr*q).ravel();z=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
    weight=(np.pi*mesh.dr*r.reshape(mesh.nr,-1)*w).ravel()[:,None]*np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()[None]
    nx,dx,ex=direct_basis(.5*mesh.rf**2,.5*r*r,degree);nz,dz,ez=direct_basis(mesh.zf,z,degree,True)
    nodes=(mesh.nr+1)*mesh.nz;theta=(mesh.nr-1)*mesh.nz
    def field(value,magnetic=False):
        if magnetic:
            p=value[:nodes].reshape(mesh.nr+1,mesh.nz);gg=value[nodes:].reshape(mesh.nr,mesh.nz)
            return np.stack((-(nx@p@dz.T)/r[:,None],r[:,None]*(ex@gg@ez.T),dx@p@nz.T),axis=-1)
        er=value[:mesh.nr*mesh.nz].reshape(mesh.nr,mesh.nz)
        et=np.pad(value[mesh.nr*mesh.nz:mesh.nr*mesh.nz+theta].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
        axial=value[mesh.nr*mesh.nz+theta:].reshape(mesh.nr+1,mesh.nz)
        return np.stack((r[:,None]*(ex@er@nz.T),(nx@(mesh.rf[:,None]*et)@nz.T)/r[:,None],nx@axial@ez.T),axis=-1)
    phase=(z+.5*d.length)/d.length;freq=np.fft.fftfreq(d.nz)*d.nz
    transform=np.exp(2j*np.pi*phase[:,None]*freq[None,:])
    if d.nz%2==0:transform[:,d.nz//2]=np.cos(np.pi*d.nz*phase)
    element=np.minimum((r*d.elements/d.radius).astype(int),d.elements-1)
    def scalar(values,transverse=False):
        radial=np.empty((len(r),d.nz))
        for e in range(d.elements):
            rows=np.flatnonzero(element==e);begin=1 if transverse and e==0 else 0
            radial[rows]=root_basis(d.x[e,begin:,0],.5*r[rows]**2)@values[e,begin:]
        result=(np.fft.fft(radial,axis=1)/d.nz)@transform.T
        if np.max(abs(result.imag))>1e-11*max(np.max(abs(result.real)),np.finfo(float).tiny):raise ValueError('Nonreal independent periodic reconstruction')
        return result.real
    u=arrays['state'];n=scalar(u[...,0]);nt=scalar(arrays['fluid_rate'][...,0])
    if np.any(n<=0):raise ValueError('Nonpositive independently reconstructed density')
    v=np.stack((r[:,None]*scalar(u[...,1]/u[...,0],True),r[:,None]*scalar(u[...,2]/u[...,0],True),scalar(u[...,3])/n),axis=-1)
    j=field(arrays['current']);b=field(arrays['magnetic'],True);b[...,2]+=ref.axial_bias
    result=dict(current_energy=float(np.sum(weight*.5*ref.reduced_mass*np.sum(j*j,axis=-1)/(E**2*n))),
        density_exchange=float(np.sum(weight*.5*ref.reduced_mass*np.sum(j*j,axis=-1)*nt/(E**2*n*n))),
        lorentz_work=float(np.sum(weight*np.sum(v*np.cross(j,b),axis=-1))))
    for prefix in ('','coupled_'):
        jt=field(arrays[prefix+'current_rate']);bt=field(arrays[prefix+'magnetic_rate'],True)
        current_rate=ref.reduced_mass/(E**2*n)*(np.sum(j*jt,axis=-1)-.5*np.sum(j*j,axis=-1)*nt/n)
        magnetic_rate=np.sum(b*bt,axis=-1)/MU0
        result[prefix+'current_energy_rate']=float(np.sum(weight*current_rate));result[prefix+'magnetic_energy_rate']=float(np.sum(weight*magnetic_rate))
    return result
