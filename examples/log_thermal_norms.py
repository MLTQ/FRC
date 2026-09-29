"""Independent physical pressure, entropy/proxy rates and local work integrals."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .physical_scalar_sampler import PhysicalScalarSampler
from .log_thermal_cases import target


def measure(mesh,reference,spec,arrays,*,order=12):
    grid=PhysicalScalarSampler(mesh,spec.get('degree',5),order)
    n=grid.density(arrays['n_cell']);nt=grid.density(arrays['nt_cell'])
    metrics={};cells={};integrals={};nodal={}
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij')
    x=.5*mesh.rf**2;dual=np.r_[x[0],.5*(x[:-1]+x[1:]),x[-1]]
    volumes=np.broadcast_to(2*np.pi*np.diff(dual)[:,None]*mesh.dz,r.shape)
    for label,gamma in (('electron',5/3),('ion',1.4)):
        s=grid.scalar(arrays[label+'_entropy']);st=grid.scalar(arrays[label+'_rate'])
        proxy=n*E*np.exp(s);p=(n/1e20)**(gamma-1)*proxy;internal=p/(gamma-1)
        expected=target(reference,grid.r,grid.z,gamma,spec.get('kind','smooth'))
        values=dict(pressure=p,temperature=p/(E*n),entropy_rate=st,proxy_rate=proxy*st,
            electric=grid.electric(arrays[label+'_electric']),internal=internal,
            internal_entropy_rate=internal*st,internal_density_rate=gamma*internal*nt/n)
        metrics[label]={k:grid.metric(value,expected[k]) for k,value in values.items()}
        metrics[label]['legacy_proxy_rate']=grid.metric(grid.scalar(arrays[label+'_legacy_rate']),expected['proxy_rate'])
        cells[label]={k:grid.cells(values[k]) for k in ('internal','internal_entropy_rate','internal_density_rate')}
        integrals[label]=dict(entropy_rate=float(np.sum(grid.weight*n*st)),
            entropy_rate_scale=float(np.sum(grid.weight*n*abs(st))),
            proxy_rate=float(np.sum(grid.weight*proxy*st)),
            proxy_rate_scale=float(np.sum(grid.weight*proxy*abs(st))))
        exact_nodes=target(reference,r,z,gamma,spec.get('kind','smooth'))['entropy_rate']
        nodal[label]={}
        for scheme,value in (('physical',arrays[label+'_rate']),
                             ('legacy',arrays[label+'_legacy_rate']/arrays[label+'_legacy_proxy'])):
            err=(value-exact_nodes)**2;norm=exact_nodes**2
            result={}
            for name,region in (('all',slice(None)),('axis_four',slice(0,4)),('interior',slice(4,-4)),('wall_four',slice(-4,None))):
                numerator=float(np.sum(volumes[region]*err[region]));denominator=float(np.sum(volumes[region]*norm[region]))
                result[name]=dict(relative_L2=float(np.sqrt(numerator/denominator)) if denominator>0 else None,
                    rms_error=float(np.sqrt(numerator/volumes[region].sum())),
                    sampled_maximum=float(np.sqrt(err[region].max())))
            nodal[label][scheme]=result
    j=grid.electric(arrays['current']);jt=grid.electric(arrays['joint_current_rate'])
    b=grid.magnetic(arrays['magnetic']);b[...,2]+=reference.axial_bias
    bt=grid.magnetic(arrays['joint_magnetic_rate'])
    coef=reference.reduced_mass/(n*E**2)
    kinetic=coef*(np.sum(j*jt,axis=-1)-.5*np.sum(j*j,axis=-1)*nt/n)
    density=-.5*coef*np.sum(j*j,axis=-1)*nt/n
    cells['field']=dict(current_energy_rate=grid.cells(kinetic),density_energy_rate=grid.cells(density),
                        magnetic_energy_rate=grid.cells(np.sum(b*bt,axis=-1)/MU0))
    metrics['joint_magnetic_rate_absolute']=grid.metric(bt,np.zeros_like(bt))
    return dict(metrics=metrics,integrals=integrals,nodal=nodal),cells
