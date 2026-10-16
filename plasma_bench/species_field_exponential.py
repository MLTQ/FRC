"""Exponential action for the frozen-population Lorentz/Ampere source."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply


def advance(momentum,electric,plasma_frequency,cyclotron,dt):
    """Evolve energy-scaled charged momenta and E without changing their norm."""
    count=len(plasma_frequency);indices=np.arange(3*count).reshape(count,3)
    field=np.arange(3*count,3*count+3);rows=[];cols=[];values=[]
    def entries(r,c,v):
        rows.extend(np.ravel(r));cols.extend(np.ravel(c));values.extend(np.ravel(v))
    for a,b,k,sign in ((0,1,2,1),(0,2,1,-1),(1,0,2,-1),(1,2,0,1),(2,0,1,1),(2,1,0,-1)):
        entries(indices[:,a],indices[:,b],sign*cyclotron[:,k])
    entries(indices,np.broadcast_to(field,indices.shape),np.repeat(plasma_frequency[:,None],3,axis=1))
    entries(np.broadcast_to(field,indices.shape),indices,-np.repeat(plasma_frequency[:,None],3,axis=1))
    operator=coo_matrix((dt*np.asarray(values),(rows,cols)),shape=(3*count+3,3*count+3)).tocsr()
    final=expm_multiply(operator,np.r_[momentum.ravel(),electric],traceA=0.)
    if not np.isfinite(final).all():raise ValueError('Nonfinite exponential electromagnetic update')
    return final[:-3].reshape(count,3),final[-3:]
