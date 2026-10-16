"""Common conservative density-cut selection and observed Gauss-field error."""
import numpy as np
from scipy.constants import elementary_charge as EV,epsilon_0 as EPS0


def select_partition(edges,present,*,roundoff_ulps=0,minimum_width=0.):
    """Remove close internal cuts; every species support endpoint is protected."""
    if not isinstance(roundoff_ulps,int) or roundoff_ulps<0:raise ValueError('Nonnegative integer ULP tolerance required')
    if not np.isfinite(minimum_width) or minimum_width<0:raise ValueError('Nonnegative finite density width required')
    protected=np.r_[True,np.any(present[:-1]!=present[1:],axis=1),True]
    kept=list(range(len(edges)))
    def remove(tolerance):
        i=0
        while tolerance and i<len(kept)-1:
            a,b=kept[i:i+2]
            if edges[b]-edges[a]<=tolerance:
                if not protected[b]:kept.pop(i+1);continue
                if not protected[a]:kept.pop(i);i=max(0,i-1);continue
            i+=1
    tolerance=roundoff_ulps*np.spacing(np.max(abs(edges)));remove(tolerance);after_roundoff=len(kept)
    remove(minimum_width)
    return kept,dict(removed_roundoff_edges=len(edges)-after_roundoff,removed_density_edges=after_roundoff-len(kept),
        geometric_tolerance_m=float(tolerance),minimum_density_width_m=float(minimum_width))


def field_projection_error(edges,amounts,charges,left_electric,model,state):
    """Observe old/new represented fields; never compensate their energy change."""
    old=np.r_[left_electric,left_electric+np.cumsum(EV*(amounts@charges))/EPS0]
    new_edges,new,_,new_energy=model.field(state)
    at_old=np.interp(np.asarray(edges,float),np.asarray(new_edges,float),np.asarray(new,float));difference=at_old-old;h=np.diff(edges)
    old_energy=EPS0/6*np.sum(h*(old[:-1]**2+old[:-1]*old[1:]+old[1:]**2))
    error2=np.sum(h*(difference[:-1]**2+difference[:-1]*difference[1:]+difference[1:]**2))/3
    return dict(projection_field_energy_change_j_m2=float(new_energy-old_energy),
        projection_field_maximum_change_v_m=float(abs(difference).max()),projection_field_l2_change_v_sqrt_m=float(np.sqrt(error2)))
