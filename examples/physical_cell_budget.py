"""Original physical cell balances and explicit one-sided face jumps."""
import numpy as np
from .current_transport_sampling import TransportSampler
from .physical_trace_sampling import TraceSampler
from .physical_parent_flux import fluxes,rates
from .matched_mode_rates import fingerprint


NAMES=('particles','current','electron_entropy','ion_entropy','energy')


def measure(d,system,state,field,rate,order):
    mesh=system.geometry.mesh;sample=TransportSampler(d,mesh,order)
    from plasma_bench.regular_current import RegularCurrent
    model=RegularCurrent(system.mapping,state[...,0],system.ion_mass)
    magnetic=field.packed_increment(system.geometry);current=model.current(magnetic)
    conserved,power,divj,gradient_square=rates(sample,state,magnetic,current,rate,field.axial_bias,system.ion_mass)
    integrate=lambda value: np.sum((sample.weight[...,None]*value).reshape(mesh.nr,order,mesh.nz,order,-1),axis=(1,3))
    volume=sample.weight.reshape(mesh.nr,order,mesh.nz,order).sum(axis=(1,3))
    volume_rate=integrate(conserved)
    q,w=np.polynomial.legendre.leggauss(order);rc=np.repeat(np.arange(mesh.nr),order);zc=np.repeat(np.arange(mesh.nz),order)
    radial=[];axial=[];raw={};fingerprints={};radial_jump=None;axial_jump=None
    wr=(np.pi*mesh.dr*sample.r.reshape(mesh.nr,order)*w).ravel()
    wz=np.broadcast_to(.5*mesh.dz*w,(mesh.nz,order)).ravel()
    for side,offset in [('left',-1),('right',0)]:
        rcell=np.minimum(np.arange(1,mesh.nr+1)+offset,mesh.nr-1)
        sr=TraceSampler(d,mesh,mesh.rf[1:],sample.z,rcell,zc)
        fr=fluxes(sr,state,magnetic,rate.electric,current,field.axial_bias,system.ion_mass)[...,0]
        rweight=2*np.pi*mesh.rf[1:,None]*wz[None,:]
        integrated=np.sum((rweight[...,None]*fr).reshape(mesh.nr,mesh.nz,order,len(NAMES)),axis=2)
        radial.append(np.concatenate((np.zeros((1,mesh.nz,len(NAMES))),integrated),axis=0))
        sz=TraceSampler(d,mesh,sample.r,mesh.zf[:-1],rc,np.arange(mesh.nz)+offset)
        fz=fluxes(sz,state,magnetic,rate.electric,current,field.axial_bias,system.ion_mass)[...,2]
        axial.append(np.sum((wr[:,None,None]*fz).reshape(mesh.nr,order,mesh.nz,len(NAMES)),axis=1))
        fingerprints[side+'_radial_flux']=fingerprint(fr);fingerprints[side+'_axial_flux']=fingerprint(fz)
        if side=='left':radial_jump=fr;axial_jump=fz
        else:radial_jump=radial_jump-fr;axial_jump=axial_jump-fz
    out=radial[0][1:]-radial[1][:-1]+np.roll(axial[0],-1,axis=1)-axial[1]
    residual=volume_rate+out
    wall=radial[0][-1].sum(axis=0)
    jump=(radial[0][1:-1]-radial[1][1:-1]).sum(axis=(0,1))+(axial[0]-axial[1]).sum(axis=(0,1))
    jump_l1=np.sum(rweight[:-1,:,None]*abs(radial_jump[:-1]),axis=(0,1))+np.sum(wr[:,None,None]*abs(axial_jump),axis=(0,1))
    closure=out.sum(axis=(0,1))-wall-jump
    charge_scale=float(np.sum(sample.weight*gradient_square))
    summaries={}
    for k,name in enumerate(NAMES):
        numerator=float(np.sum(residual[...,k]**2/volume))
        denominator=charge_scale if name=='current' else float(np.sum((volume_rate[...,k]**2+out[...,k]**2)/volume))
        summaries[name]=dict(cell_balance_relative_L2=float(np.sqrt(numerator/max(denominator,np.finfo(float).tiny))),
            cell_balance_absolute_L2=float(np.sqrt(numerator)),cell_balance_max_density=float(np.max(abs(residual[...,k]/volume))),
            integrated_time_rate=float(volume_rate[...,k].sum()),integrated_cell_outflux=float(out[...,k].sum()),
            exterior_wall_flux=float(wall[k]),internal_outward_jump=float(jump[k]),
            internal_trace_jump_L1=float(jump_l1[k]),geometric_flux_closure=float(closure[k]),
            global_rate_plus_exterior=float(volume_rate[...,k].sum()+wall[k]))
    integrated_divj=integrate(divj[...,None])[...,0]
    powers={name:float(np.sum(sample.weight*value)) for name,value in power.items()}
    raw.update(volume=volume,volume_rate=volume_rate,cell_outflux=out,cell_balance=residual,
        radial_left_integrals=radial[0],radial_right_integrals=radial[1],
        axial_left_integrals=axial[0],axial_right_integrals=axial[1],strong_current_divergence_integrals=integrated_divj)
    return dict(order=order,channels=summaries,continuous_power_components=powers,
        strong_current_divergence_relative_L2=float(np.sqrt(np.sum(sample.weight*divj*divj)/max(charge_scale,np.finfo(float).tiny))),
        current_cell_divergence_theorem_max=float(np.max(abs(integrated_divj-out[...,1]))),
        physical_fingerprints=fingerprints,
        axis_flux='zero cylindrical area',axial_boundary='periodic seam included among internal jumps',
        current_normalization='physical gradient Frobenius L2; no evolved charge state'),raw
