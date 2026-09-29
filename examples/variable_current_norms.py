"""Independent primitive-polynomial density and physical field/energy audit."""
import numpy as np
from numpy.polynomial import Polynomial
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .compatible_field_norms import direct_basis


def primitive_map(edges,points,degree,periodic=False):
    """Differentiate root-product interpolants of exact cumulative cell integrals."""
    count=len(edges)-1;out=np.zeros((len(points),count))
    for row,point in enumerate(points):
        cell=min(np.searchsorted(edges,point,side='right')-1,count-1)
        for offset in ((degree+1)//2,(degree-1)//2):
            start=cell-offset
            if not periodic:start=min(max(start,0),count-degree-1)
            ids=np.arange(start,start+degree+2)
            coords=edges[0]+ids*(edges[1]-edges[0]) if periodic else edges[ids]
            scale=coords[-1]-coords[0];roots=(coords-point)/scale
            derivative=[]
            for i in range(len(ids)):
                others=np.delete(roots,i)
                poly=Polynomial.fromroots(others)/np.prod(roots[i]-others)
                derivative.append(poly.deriv()(0.)/scale)
            for i,index in enumerate(ids[:-1]):
                out[row,index%count]+=.5*(coords[i+1]-coords[i])*sum(derivative[i+1:])
    return out


def measure(mesh,reference,degree,arrays,*,order=12):
    q,w=np.polynomial.legendre.leggauss(order)
    r=(mesh.r[:,None]+.5*mesh.dr*q).ravel();z=(mesh.z[:,None]+.5*mesh.dz*q).ravel()
    nx,dx,ex=direct_basis(.5*mesh.rf**2,.5*r*r,degree)
    nz,dz,ez=direct_basis(mesh.zf,z,degree,True)
    rx=primitive_map(.5*mesh.rf**2,.5*r*r,degree)
    rz=primitive_map(mesh.zf,z,degree,True)
    weight=(np.pi*mesh.dr*r.reshape(mesh.nr,-1)*w).ravel()[:,None]*np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()[None]
    nodes=(mesh.nr+1)*mesh.nz;theta=(mesh.nr-1)*mesh.nz
    def field(value,magnetic=False):
        if magnetic:
            p=value[:nodes].reshape(mesh.nr+1,mesh.nz);g=value[nodes:].reshape(mesh.nr,mesh.nz)
            return np.stack((-(nx@p@dz.T)/r[:,None],r[:,None]*(ex@g@ez.T),dx@p@nz.T),axis=-1)
        er=value[:mesh.nr*mesh.nz].reshape(mesh.nr,mesh.nz)
        et=np.pad(value[mesh.nr*mesh.nz:mesh.nr*mesh.nz+theta].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
        axial=value[mesh.nr*mesh.nz+theta:].reshape(mesh.nr+1,mesh.nz)
        return np.stack((r[:,None]*(ex@er@nz.T),(nx@(mesh.rf[:,None]*et)@nz.T)/r[:,None],nx@axial@ez.T),axis=-1)
    exact=reference.fields(r[:,None],z[None]);metrics={}
    values={key:field(arrays[key],key.startswith('magnetic')) for key in ('magnetic','magnetic_rate','current','current_rate','electric')}
    values['n']=rx@arrays['n_cell']@rz.T;values['nt']=rx@arrays['nt_cell']@rz.T
    for label,gamma in (('electron',5/3),('ion',1.4)):
        values[label+'_electric']=field(arrays[label+'_electric'])
        exact[label+'_electric']=reference.thermal(r[:,None],z[None],gamma=gamma)['electric']
    n=values['n'];nt=values['nt'];j=values['current'];jt=values['current_rate']
    coefficient=reference.reduced_mass/(n*E**2)
    full_b=values['magnetic'].copy();full_b[...,2]+=reference.axial_bias
    values['current_energy']=.5*coefficient*np.sum(j*j,axis=-1)
    values['density_energy_rate']=-.5*coefficient*np.sum(j*j,axis=-1)*nt/n
    values['current_energy_rate']=coefficient*np.sum(j*jt,axis=-1)+values['density_energy_rate']
    values['magnetic_energy_rate']=np.sum(full_b*values['magnetic_rate'],axis=-1)/MU0
    for key,value in values.items():
        wanted=exact[key];error=(value-wanted)**2;norm=wanted**2
        if error.ndim==3:error=error.sum(axis=-1);norm=norm.sum(axis=-1)
        metrics[key]=dict(relative_L2=float(np.sqrt(np.sum(weight*error)/np.sum(weight*norm))),
            rms_error=float(np.sqrt(np.sum(weight*error)/mesh.volumes.sum())),
            axis_relative_L2=float(np.sqrt(np.sum(weight[:order]*error[:order])/np.sum(weight[:order]*norm[:order]))),
            sampled_maximum=float(np.sqrt(error.max())))
    def cells(value):
        return np.sum((weight*value).reshape(mesh.nr,order,mesh.nz,order),axis=(1,3))/mesh.volumes
    energies={key:cells(values[key]) for key in ('current_energy','current_energy_rate','density_energy_rate','magnetic_energy_rate')}
    return dict(metrics=metrics,integrals={k:float(np.sum(mesh.volumes*v)) for k,v in energies.items()}),energies
