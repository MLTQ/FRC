"""Independent physical current and one-sided magnetic/current face errors."""
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.regular_current import RegularCurrent
from .current_transport_sampling import TransportSampler
from .physical_trace_sampling import TraceSampler
from .contact_mode_fields import magnetic_curl
from .regular_plasma_norms import norm
from .matched_mode_rates import fingerprint


def measure(mode,d,system,state,magnetic,order):
    mesh=system.geometry.mesh;sample=TransportSampler(d,mesh,order)
    model=RegularCurrent(system.mapping,state[...,0],system.ion_mass);coefficients=model.current(magnetic)
    current,gradient=sample.electric(coefficients,True);target=mode.fields(sample.R,sample.Z)
    strong=magnetic_curl(sample,magnetic)/MU0;divergence=np.trace(gradient,axis1=-2,axis2=-1)
    volume=dict(initial_current=norm(sample,current,target['current']),
        strong_ampere=norm(sample,current,strong),strong_current_target=norm(sample,strong,target['current']),
        current_divergence_relative=float(np.sqrt(np.sum(sample.weight*divergence**2)/np.sum(sample.weight*np.sum(gradient**2,axis=(-2,-1))))))
    q,w=np.polynomial.legendre.leggauss(order);rc=np.repeat(np.arange(mesh.nr),order);zc=np.repeat(np.arange(mesh.nz),order)
    wr=(np.pi*mesh.dr*sample.r.reshape(mesh.nr,order)*w).ravel()
    wz=np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()
    summaries={};raw={};prints={}
    for axis in ('radial','axial'):
        normal=0 if axis=='radial' else 2;tangential=[1,2] if axis=='radial' else [0,1]
        traces=[]
        for offset in (-1,0):
            if axis=='radial':
                # Internal faces only; the closed exterior wall is a separate trace.
                s=TraceSampler(d,mesh,mesh.rf[1:-1],sample.z,np.arange(1,mesh.nr)+offset,zc)
                weight=2*np.pi*mesh.rf[1:-1,None]*wz[None,:]
            else:
                s=TraceSampler(d,mesh,sample.r,mesh.zf[:-1],rc,np.arange(mesh.nz)+offset)
                weight=np.broadcast_to(wr[:,None],(len(sample.r),mesh.nz))
            analytic=mode.fields(s.R,s.Z)
            values=dict(normal_current=s.electric(coefficients)[...,normal:normal+1],tangential_B=s.magnetic(magnetic)[...,tangential])
            expected=dict(normal_current=analytic['current'][...,normal:normal+1],tangential_B=analytic['physical'][...,6:9][...,tangential])
            traces.append(values)
            for name,value in values.items():prints[f'{axis}_{offset}_{name}']=fingerprint(value)
        metrics={}
        for name in traces[0]:
            denominator=np.sum(weight[...,None]*expected[name]**2)
            jump=traces[0][name]-traces[1][name]
            metrics[name]=dict(left_relative_L2=float(np.sqrt(np.sum(weight[...,None]*(traces[0][name]-expected[name])**2)/denominator)),
                right_relative_L2=float(np.sqrt(np.sum(weight[...,None]*(traces[1][name]-expected[name])**2)/denominator)),
                jump_over_reference_L2=float(np.sqrt(np.sum(weight[...,None]*jump**2)/denominator)),
                jump_max=float(np.max(abs(jump))))
            raw[f'{axis}_{name}_jump_squared_integrals']=np.sum(weight[...,None]*jump**2,axis=1 if axis=='radial' else 0)
        summaries[axis]=metrics
    return dict(order=order,volume=volume,internal_faces=summaries,physical_fingerprints=prints),raw
