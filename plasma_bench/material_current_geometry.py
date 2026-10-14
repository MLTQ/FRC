"""Physical meridional material metrics shared by isentropic and thermal current."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse import diags
from scipy.sparse.linalg import splu


def meridional_geometry(model,state,fields=None):
    """No thermal restriction; still the toroidal-field, nonrotating sector."""
    f=model.fluid;space=model.space;g=f.fields(state) if fields is None else fields
    if np.any(state.angle!=0) or np.any(state.angular_momentum!=0):
        raise ValueError('Toroidal-field invariant sector requires zero bulk rotation')
    deformation=np.einsum('cqnb,cna->cqab',f.reference_geometry['gradient'],state.coordinates[f.ids]).reshape(-1,2,2)
    inverse=np.linalg.inv(deformation);volume=g['volume'].ravel();r=g['points'][...,0].ravel();ratio=volume/model.volume0
    er=diags(inverse[:,0,0])@space.er+diags(inverse[:,1,0])@space.ez
    ez=diags(inverse[:,0,1])@space.er+diags(inverse[:,1,1])@space.ez
    factor=r/space.reference_r/ratio
    bmap=diags(factor)@space.magnetic;curl=diags(factor)@space.curl
    mass=(er.T@diags(model.particles)@er+ez.T@diags(model.particles)@ez).tocsc();solve=splu(mass).solve
    h=(bmap.T@diags(volume/MU0)@bmap).toarray();d=-(bmap.T@diags(volume/MU0)@curl).tocsr()
    gradient=np.einsum('cqnb,cna->cqab',g['gradient'],g['velocity'][f.ids]).reshape(-1,2,2)
    radial_rate=g['sampled_velocity'][...,0].ravel()/r;divergence=np.trace(gradient,axis1=1,axis2=2)+radial_rate
    metric_rate=2*radial_rate-divergence
    hdot=(bmap.T@diags(volume*metric_rate/MU0)@bmap).toarray()
    ddot=-(bmap.T@diags(volume*metric_rate/MU0)@curl).tocsr()
    grdot=-diags(gradient[:,0,0])@er-diags(gradient[:,1,0])@ez
    gzdot=-diags(gradient[:,0,1])@er-diags(gradient[:,1,1])@ez
    return dict(fields=g,er=er,ez=ez,bmap=bmap,curl=curl,mass=mass,solve=solve,h=h,d=d,hdot=hdot,ddot=ddot,
        grdot=grdot,gzdot=gzdot,gradient=gradient,volume=volume,r=r)
