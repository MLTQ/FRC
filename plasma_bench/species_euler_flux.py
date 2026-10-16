"""Conservative HLLC flux for species moments in a fixed velocity frame."""
import numpy as np


def physical(primitive, mass_ratio, frame_normal, axis, gamma):
    number=primitive[...,0];pressure=primitive[...,1];velocity=primitive[...,2:]
    if np.any(number<0) or np.any(pressure<0) or np.any((number==0)&(pressure!=0)):
        raise ValueError('Physical face density/pressure required; no floor')
    q=np.empty_like(primitive);q[...,0]=number
    q[...,1:4]=(mass_ratio[...,None]*velocity)*number[...,None]
    q[...,4]=pressure/(gamma-1)+(.5*mass_ratio*np.sum(velocity**2,axis=-1))*number
    normal=velocity[...,axis]+frame_normal
    flux=normal[...,None]*q;flux[...,1+axis]+=pressure;flux[...,4]+=pressure*velocity[...,axis]
    sound=np.sqrt(gamma*np.divide(pressure,number,out=np.zeros_like(number),where=number>0)/mass_ratio)
    return q,flux,normal,sound


def hllc(left, right, mass_ratio, frame_normal, *, axis=2, gamma=5/3):
    """Dimensionless face flux; both exact empty sides return exactly zero.

    Primitive slots are number, pressure, and three relative velocities. HLL
    handles a vacuum side or an inadmissible HLLC star state, without state repair.
    """
    ql,fl,ul,cl=physical(left,mass_ratio,frame_normal,axis,gamma)
    qr,fr,ur,cr=physical(right,mass_ratio,frame_normal,axis,gamma)
    nl,nr=left[...,0],right[...,0];root_l,root_r=np.sqrt(nl),np.sqrt(nr)
    denominator=root_l+root_r
    weight=np.divide(root_l,denominator,out=np.zeros_like(nl),where=denominator>0)
    roe_velocity=weight*ul+(1-weight)*ur
    roe_sound=np.sqrt(weight*cl**2+(1-weight)*cr**2+
        .5*(gamma-1)*weight*(1-weight)*np.sum((left[...,2:]-right[...,2:])**2,axis=-1))
    sl=np.minimum(np.minimum(ul-cl,ur-cr),roe_velocity-roe_sound)
    sr=np.maximum(np.maximum(ul+cl,ur+cr),roe_velocity+roe_sound)
    sl=np.where((nl==0)&(nr>0),np.minimum(sl,ur-2*cr/(gamma-1)),sl)
    sr=np.where((nr==0)&(nl>0),np.maximum(sr,ul+2*cl/(gamma-1)),sr)
    result=np.where((sl>=0)[...,None],fl,fr).copy()
    fan=(sl<0)&(sr>0)
    if not np.any(fan):return result
    a,b=sl[fan],sr[fan];l,r=left[fan],right[fan];L,R=ql[fan],qr[fan]
    FL,FR=fl[fan],fr[fan];u,v=ul[fan],ur[fan]
    hll=(b[:,None]*FL-a[:,None]*FR+(a*b)[:,None]*(R-L))/(b-a)[:,None]
    scale=np.maximum(l[:,0],r[:,0])
    masses=np.broadcast_to(mass_ratio,nl.shape)[fan]
    rl=(l[:,0]/scale)*masses;rr=(r[:,0]/scale)*masses
    pl,pr=l[:,1]/scale,r[:,1]/scale
    contact=(pr-pl+rl*u*(a-u)-rr*v*(b-v))/(rl*(a-u)-rr*(b-v))
    pressure=(pl+rl*(a-u)*(contact-u))*scale
    choose=contact>=0
    Q=np.where(choose[:,None],L,R);F=np.where(choose[:,None],FL,FR)
    primitive=np.where(choose[:,None],l,r);s=np.where(choose,a,b);normal=np.where(choose,u,v)
    valid=(l[:,0]>0)&(r[:,0]>0)&(pressure>=0)&np.isfinite(contact)&(s!=contact)
    indices=np.flatnonzero(valid)
    if len(indices):
        k=indices;factor=(s[k]-normal[k])/(s[k]-contact[k])
        star=Q[k].copy();star[:,0]=primitive[k,0]*factor
        velocity=primitive[k,2:].copy();velocity[:,axis]=contact[k]-frame_normal
        star[:,1:4]=(masses[k,None]*velocity)*star[:,0,None]
        star[:,4]=((s[k]-normal[k])*Q[k,4]-primitive[k,1]*primitive[k,2+axis]
            +pressure[k]*(contact[k]-frame_normal))/(s[k]-contact[k])
        internal=star[:,4]-(.5*masses[k]*np.sum(velocity**2,axis=1))*star[:,0]
        admissible=(star[:,0]>0)&(internal>=0)&np.isfinite(star).all(axis=1)
        kk=k[admissible]
        hll[kk]=F[kk]+s[kk,None]*(star[admissible]-Q[kk])
    result[fan]=hll
    if not np.isfinite(result).all():raise ValueError('Nonfinite species flux; no repair')
    return result
