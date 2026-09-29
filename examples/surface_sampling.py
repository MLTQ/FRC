"""Evaluate saved Q1 material states at common physical reference coordinates."""
import numpy as np


def locate(arrays,configuration,points):
    """Invert the initial straight radial-sector map, including graded shells."""
    nt=configuration["angular_cells"];nr=configuration["radial_cells"]
    x=np.asarray(points,float).reshape(-1,2);theta=np.linspace(0.,np.pi,nt+1)
    angle=np.arctan2(x[:,0],-x[:,1]);j=np.clip(np.searchsorted(theta,angle,side="right")-1,0,nt-1)
    e0=np.column_stack((np.sin(theta[j]),-np.cos(theta[j])))
    e1=np.column_stack((np.sin(theta[j+1]),-np.cos(theta[j+1])))
    det=np.sin(np.pi/nt)
    c0=(x[:,0]*e1[:,1]-e1[:,0]*x[:,1])/det
    c1=(e0[:,0]*x[:,1]-x[:,0]*e0[:,1])/det
    radius=c0+c1
    shells=np.r_[0.,np.linalg.norm(arrays["t0_reference_coordinates"][1::nt+1],axis=1)]
    if len(shells)!=nr+1 or np.any(radius < -2e-12) or np.any(radius>shells[-1]+2e-12):
        raise ValueError("Comparison point is outside the reference material polygon")
    radius=np.clip(radius,0.,shells[-1])
    k=np.clip(np.searchsorted(shells,radius,side="right")-1,0,nr-1)
    xi=(radius-shells[k])/(shells[k+1]-shells[k])
    tau=np.divide(c1,radius,out=np.zeros_like(radius),where=radius>0)
    tau=np.clip(tau,0.,1.)
    shape=np.stack(((1-xi)*(1-tau),xi*(1-tau),xi*tau,(1-xi)*tau),axis=-1)
    derivative=np.stack((np.stack((-(1-tau),1-tau,tau,-tau),axis=-1),
                         np.stack((-(1-xi),-xi,xi,1-xi),axis=-1)),axis=-1)
    cells=k*nt+j;ids=arrays["t0_cells"][cells]
    return cells,ids,shape,derivative


def nodal(arrays,configuration,prefix,points,key):
    _,ids,shape,_=locate(arrays,configuration,points)
    return np.einsum("qn,qn...->q...",shape,arrays[prefix+"_"+key][ids])


def fields(arrays,configuration,prefix,points):
    cells,ids,shape,deriv=locate(arrays,configuration,points)
    x=arrays[prefix+"_coordinates"][ids]
    position=np.einsum("qn,qna->qa",shape,x)
    jac=np.einsum("qnb,qna->qab",deriv,x)
    initial_jac=np.einsum("qnb,qna->qab",deriv,arrays["t0_reference_coordinates"][ids])
    gradient=np.einsum("qnb,qba->qna",deriv,np.linalg.inv(jac))
    a=arrays[prefix+"_total_material_a"][ids]
    value=np.einsum("qn,qn->q",shape,a);da=np.einsum("qna,qn->qa",gradient,a)
    r=position[:,0]
    return {"position":position,"volume_ratio":r/np.asarray(points)[:,0]*np.linalg.det(jac)/np.linalg.det(initial_jac),
            "velocity":np.einsum("qn,qna->qa",shape,arrays[prefix+"_velocity"][ids]),
            "magnetic":np.column_stack((-r*da[:,1],2*value+r*da[:,0])),
            "current":r*np.einsum("qn,qn->q",shape,arrays[prefix+"_current_over_radius"][ids]),
            "density":arrays[prefix+"_density"][cells],"pressure":arrays[prefix+"_pressure"][cells],
            "heat_density":(arrays[prefix+"_cell_ohmic_heat"]/arrays[prefix+"_volume"])[cells]}
