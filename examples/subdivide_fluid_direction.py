"""Nested local subdivision of an archived conservative fluid direction."""
import numpy as np


def transfer(source,target,value):
    if target.elements%source.elements or target.degree!=source.degree or target.nz!=source.nz:
        raise ValueError('Nested radial elements and unchanged degree/axial grid required')
    if not np.array_equal(target.z,source.z) or target.radius!=source.radius:
        raise ValueError('Original physical domain and axial nodes required')
    if np.shape(value)!=source.shape:raise ValueError('Complete conservative direction required')
    ratio=target.elements//source.elements;result=np.zeros(target.shape)
    for child in range(target.elements):
        parent=child//ratio
        for slot in range(6):
            begin=1 if parent==0 and slot in (1,2) else 0
            nodes=source.x[parent,begin:,0];scale=nodes[-1]-nodes[0]
            x=(nodes-nodes[0])/scale
            for local in range(target.degree+1):
                if child==0 and local==0 and slot in (1,2):continue
                if local==0 and child%ratio==0 and begin==0:
                    result[child,local,:,slot]=value[parent,0,:,slot];continue
                if local==target.degree and child%ratio==ratio-1:
                    result[child,local,:,slot]=value[parent,-1,:,slot];continue
                point=(target.x[child,local,0]-nodes[0])/scale
                weights=np.ones(len(nodes))
                for i in range(len(nodes)):
                    for j in range(len(nodes)):
                        if i!=j:weights[i]*=(point-x[j])/(x[i]-x[j])
                result[child,local,:,slot]=weights@value[parent,begin:,:,slot]
    if not np.isfinite(result).all():raise ValueError('Finite transferred direction required')
    return result


def parent_faces(value,parents):
    if value.shape[0]%parents:raise ValueError('Nested parent count required')
    ratio=value.shape[0]//parents
    return np.stack((value[::ratio,0],value[ratio-1::ratio,-1]),axis=1)
