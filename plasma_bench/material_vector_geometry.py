"""Full axisymmetric Piola metrics, including material rotation and winding."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse import diags
from scipy.sparse.linalg import splu


def vector_geometry(model,state):
    f=model.fluid;space=model.space;fields=f.fields(state);shape=f.particles.shape
    reference=f.reference_geometry;r=fields['points'][...,0].ravel();r0=reference['points'][...,0].ravel()
    deformation=np.zeros((len(r),3,3))
    poloidal=np.einsum('cqnb,cna->cqab',reference['gradient'],state.coordinates[f.ids]).reshape(-1,2,2)
    twist=np.einsum('cqna,cn->cqa',reference['gradient'],state.angle[f.ids]).reshape(-1,2)
    deformation[:,0,0]=poloidal[:,0,0];deformation[:,0,2]=poloidal[:,0,1]
    deformation[:,2,0]=poloidal[:,1,0];deformation[:,2,2]=poloidal[:,1,1]
    deformation[:,1,0]=r*twist[:,0];deformation[:,1,1]=r/r0;deformation[:,1,2]=r*twist[:,1]
    volume=fields['volume'].ravel();jacobian=volume/model.volume0;inverse=np.linalg.inv(deformation)
    electric=tuple(sum((diags(inverse[:,j,i])@space.electric[j] for j in range(3))) for i in range(3))
    magnetic=tuple(sum((diags(deformation[:,i,j]/jacobian)@space.magnetic[j] for j in range(3))) for i in range(3))
    curl=tuple(sum((diags(deformation[:,i,j]/jacobian)@space.curl[j] for j in range(3))) for i in range(3))
    mass=sum(g.T@diags(model.particles)@g for g in electric).tocsc();solve=splu(mass).solve
    h=sum(g.T@diags(volume/MU0)@g for g in magnetic).toarray()
    d=-sum(b.T@diags(volume/MU0)@c for b,c in zip(magnetic,curl)).tocsr()
    grad=np.einsum('cqnb,cna->cqab',fields['gradient'],fields['velocity'][f.ids]).reshape(-1,2,2)
    omega=fields['sampled_omega'].ravel()
    omega_gradient=np.einsum('cqna,cn->cqa',fields['gradient'],fields['omega'][f.ids]).reshape(-1,2)
    # The current cylindrical basis rotates with material angular speed omega.
    # Fdot_local=(grad(v)-omega*J_rotation)*F. Its symmetric part is grad(v)'s.
    rate=np.zeros_like(deformation)
    rate[:,0,0]=grad[:,0,0];rate[:,0,2]=grad[:,0,1]
    rate[:,2,0]=grad[:,1,0];rate[:,2,2]=grad[:,1,1]
    rate[:,1,0]=r*omega_gradient[:,0];rate[:,1,1]=fields['sampled_velocity'][...,0].ravel()/r;rate[:,1,2]=r*omega_gradient[:,1]
    divergence=np.trace(rate,axis1=1,axis2=2)
    metric=rate+rate.transpose(0,2,1)-divergence[:,None,None]*np.eye(3)
    hdot=sum(magnetic[i].T@diags(volume*metric[:,i,j]/MU0)@magnetic[j] for i in range(3) for j in range(3)).toarray()
    ddot=-sum(magnetic[i].T@diags(volume*metric[:,i,j]/MU0)@curl[j] for i in range(3) for j in range(3)).tocsr()
    electric_rate=tuple(-sum(diags(rate[:,j,i])@electric[j] for j in range(3)) for i in range(3))
    vorticity=np.column_stack((-r*omega_gradient[:,1],grad[:,0,1]-grad[:,1,0],2*omega+r*omega_gradient[:,0]))
    velocity=np.column_stack((fields['sampled_velocity'][...,0].ravel(),r*omega,fields['sampled_velocity'][...,1].ravel()))
    return dict(fields=fields,electric=electric,magnetic=magnetic,curl=curl,mass=mass,solve=solve,h=h,d=d,
        hdot=hdot,ddot=ddot,electric_rate=electric_rate,vorticity=vorticity,volume=volume,r=r,velocity=velocity,deformation=deformation)
