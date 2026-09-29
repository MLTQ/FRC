"""Independent direct-polynomial physical field norms for compatible cochains."""
import numpy as np
from numpy.polynomial import Polynomial
from scipy.constants import mu_0 as MU0,elementary_charge as E


def direct_basis(edges,points,degree,periodic=False):
    """Root-product Lagrange polynomials, independent of production moment solves."""
    count=len(edges)-1;size=count if periodic else count+1
    N=np.zeros((len(points),size));dN=N.copy();edge=np.zeros((len(points),count))
    for row,point in enumerate(points):
        cell=min(np.searchsorted(edges,point,side='right')-1,count-1)
        start=cell-(degree-1)//2
        if not periodic:start=min(max(start,0),count-degree)
        ids=np.arange(start,start+degree+1)
        coords=edges[0]+ids*(edges[1]-edges[0]) if periodic else edges[ids]
        scale=coords[-1]-coords[0];roots=(coords-point)/scale
        derivative=[]
        for local,index in enumerate(ids):
            others=np.delete(roots,local)
            poly=Polynomial.fromroots(others)/np.prod(roots[local]-others)
            column=index%size if periodic else index
            N[row,column]+=poly(0.);dN[row,column]+=poly.deriv()(0.)/scale
            derivative.append(poly.deriv()(0.)/scale)
        for k,index in enumerate(ids[:-1]):edge[row,index%count if periodic else index]=-sum(derivative[:k+1])*(coords[k+1]-coords[k])
    return N,dN,edge


def measure(mesh,ref,degree,arrays,*,bias_in_state=True,order=12):
    q,w=np.polynomial.legendre.leggauss(order)
    r=(mesh.r[:,None]+.5*mesh.dr*q).ravel();z=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
    nx,dx,ex=direct_basis(.5*mesh.rf**2,.5*r*r,degree)
    nz,dz,ez=direct_basis(mesh.zf,z,degree,True)
    weight=(np.pi*mesh.dr*r.reshape(mesh.nr,-1)*w).ravel()[:,None]*np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()[None]
    nodes=(mesh.nr+1)*mesh.nz;ct=(mesh.nr-1)*mesh.nz
    def magnetic(packed):
        psi=packed[:nodes].reshape(mesh.nr+1,mesh.nz);g=packed[nodes:].reshape(mesh.nr,mesh.nz)
        return np.stack((-(nx@psi@dz.T)/r[:,None],r[:,None]*(ex@g@ez.T),dx@psi@nz.T),axis=-1)
    def electric(packed):
        er=packed[:mesh.nr*mesh.nz].reshape(mesh.nr,mesh.nz)
        et=np.pad(packed[mesh.nr*mesh.nz:mesh.nr*mesh.nz+ct].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
        axial=packed[mesh.nr*mesh.nz+ct:].reshape(mesh.nr+1,mesh.nz)
        return np.stack((r[:,None]*(ex@er@nz.T),(nx@(mesh.rf[:,None]*et)@nz.T)/r[:,None],nx@axial@ez.T),axis=-1)
    exact=ref.fields(r[:,None],z[None]);target=ref.derivatives(r[:,None],z[None])[2]
    b=magnetic(arrays['magnetic']);bt=magnetic(arrays['magnetic_rate'])
    j=electric(arrays['current']);jt=electric(arrays['current_rate']);field=electric(arrays['electric_full'])
    full_b=b.copy();full_b[...,2]+=float(arrays.get('axial_bias_t',0.))
    energy=float(np.sum(weight*np.sum(full_b*full_b,axis=-1))/(2*MU0))
    kinetic=float(np.sum(weight*np.sum(j*j,axis=-1))*ref.reduced_mass/(2*ref.number_density*E**2))
    if bias_in_state and 'axial_bias_t' not in arrays:b[...,2]-=ref.bias
    metrics={}
    for key,value,wanted in (('magnetic',b,exact['magnetic']),('current',j,exact['current']),
        ('magnetic_rate',bt,target['magnetic'].real),('current_rate',jt,target['current'].real),('electric',field,exact['electric'])):
        error=np.sum((value-wanted)**2,axis=-1);norm=np.sum(wanted*wanted,axis=-1)
        metrics[key]=dict(relative_L2=float(np.sqrt(np.sum(weight*error)/np.sum(weight*norm))),
            axis_relative_L2=float(np.sqrt(np.sum(weight[:order]*error[:order])/np.sum(weight[:order]*norm[:order]))),
            rms_error=float(np.sqrt(np.sum(weight*error)/mesh.volumes.sum())),sampled_maximum=float(np.sqrt(error.max())))
    return dict(metrics=metrics,magnetic_energy_j=energy,current_energy_j=kinetic)
