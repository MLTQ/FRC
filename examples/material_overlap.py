"""Positive quadrature on intersections of two initial convex material meshes."""
import numpy as np
from numpy.polynomial.legendre import leggauss


def cross(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def intersection(subject,clip):
    """Sutherland-Hodgman clipping; both convex polygons have positive orientation."""
    result=np.asarray(subject,float)
    for start,end in zip(clip,np.roll(clip,-1,axis=0)):
        edge=end-start
        if np.all(edge==0):continue
        if len(result)<3:return np.empty((0,2))
        output=[];previous=result[-1];old=float(cross(edge,previous-start))
        for point in result:
            new=float(cross(edge,point-start))
            if (old>=0)!=(new>=0):
                output.append(previous+old/(old-new)*(point-previous))
            if new>=0:output.append(point)
            previous,old=point,new
        result=np.asarray(output).reshape(-1,2)
    return result


def intersections(a,b):
    """Intersect initial Q1 polygons, using only bounding boxes for candidate rejection."""
    pa=a["t0_reference_coordinates"][a["t0_cells"]]
    pb=b["t0_reference_coordinates"][b["t0_cells"]]
    lower=pb.min(axis=1);upper=pb.max(axis=1);polygons=[]
    scale=max(float(np.max(abs(pa))),float(np.max(abs(pb))))
    for cell in pa:
        lo=cell.min(axis=0);hi=cell.max(axis=0)
        candidates=np.where(np.all(upper>=lo,axis=1)&np.all(lower<=hi,axis=1))[0]
        for index in candidates:
            poly=intersection(cell,pb[index])
            if len(poly)<3:continue
            shifted=poly-poly[0]
            area=.5*float(np.sum(cross(shifted,np.roll(shifted,-1,axis=0))))
            if area>1e-26*scale**2:polygons.append(poly)
    return polygons


def quadrature(polygons,order=4,chunk_size=65536):
    """Triangulate each convex overlap and integrate with a positive Duffy rule."""
    q,w=leggauss(order);u,v=np.meshgrid((q+1)/2,(q+1)/2,indexing="ij")
    u,v=u.ravel(),v.ravel();base=(w[:,None]*w[None,:]).ravel()*.25*(1-u)
    points=[];weights=[];count=0
    for poly in polygons:
        center=poly.mean(axis=0)
        for first,last in zip(poly,np.roll(poly,-1,axis=0)):
            e0,e1=first-center,last-center;det=float(cross(e0,e1))
            if det<=0:continue # Repeated collinear clipping vertices have no area.
            points.append(center+u[:,None]*e0+((1-u)*v)[:,None]*e1)
            weights.append(det*base);count+=len(u)
            if count>=chunk_size:
                yield np.concatenate(points),np.concatenate(weights)
                points=[];weights=[];count=0
    if points:yield np.concatenate(points),np.concatenate(weights)
