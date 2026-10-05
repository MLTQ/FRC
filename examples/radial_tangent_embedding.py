"""Exact element subdivision of the declared conservative tangent."""
import numpy as np


def embed(old, new, value):
    if new.elements % old.elements or new.nz != old.nz or new.degree != old.degree or new.radius != old.radius:
        raise ValueError('Same degree, radius and axial grid with integer subdivision required')
    if value.shape != old.shape: raise ValueError('Complete original tangent required')
    result=np.zeros(new.shape);ratio=new.elements//old.elements
    for child in range(new.elements):
        parent=child//ratio
        for transverse in (False,True):
            begin=int(transverse and parent==0);nodes=old.x[parent,begin:,0]
            origin=nodes[0];scale=nodes[-1]-origin
            roots=(nodes-origin)/scale;points=(new.x[child,:,0]-origin)/scale
            matrix=np.column_stack([np.prod(points[:,None]-np.delete(roots,k),axis=1)/
                np.prod(roots[k]-np.delete(roots,k)) for k in range(len(roots))])
            for channel in ((1,2) if transverse else (0,3,4,5)):
                result[child,:,:,channel]=matrix@value[parent,begin:,:,channel]
    result[0,0,:,1:3]=0.
    return result
