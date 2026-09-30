"""Independent root-product flux and physical moving-thermal norms."""
import numpy as np
from scipy.constants import elementary_charge as E
from .physical_scalar_sampler import PhysicalScalarSampler


def flux_basis(edges,points,degree,periodic=False):
    """Interpolate cumulative face values, independently of production moments."""
    count=len(edges)-1;out=np.zeros((len(points),count if periodic else count+1))
    for row,point in enumerate(points):
        cell=min(np.searchsorted(edges,point,side='right')-1,count-1)
        for offset in ((degree+1)//2,(degree-1)//2):
            start=cell-offset if periodic else min(max(cell-offset,0),count-degree-1)
            ids=np.arange(start,start+degree+2)
            coords=edges[0]+ids*(edges[1]-edges[0]) if periodic else edges[ids]
            scaled=(coords-point)/(coords[-1]-coords[0])
            for i,index in enumerate(ids):
                other=np.delete(scaled,i)
                weight=.5*np.prod(-other)/np.prod(scaled[i]-other)
                out[row,index%count if periodic else index]+=weight
    return out


def rate_norms(mesh,ref,spec,arrays,order=14):
    grid=PhysicalScalarSampler(mesh,spec.get('degree',5),order);degree=spec.get('degree',5)
    hx=flux_basis(.5*mesh.rf**2,.5*grid.r[:,0]**2,degree)
    hz=flux_basis(mesh.zf,grid.z[0],degree,True)
    particle=np.stack(((hx@arrays['radial']@grid.Rz.T)/grid.r,np.zeros(np.broadcast_shapes(grid.r.shape,grid.z.shape)),
                       grid.Rx@arrays['axial']@hz.T),axis=-1)
    n=grid.density(arrays['density']);nt=grid.density(arrays['density_rate']);f=ref.fields(grid.r,grid.z)
    metrics=dict(density=grid.metric(n,f['n']),density_rate=grid.metric(nt,f['nt']),particle_flux=grid.metric(particle,f['particle']))
    local={};balances={}
    for label,gamma in (('electron',5/3),('ion',1.4)):
        s=grid.scalar(arrays[label]);st=grid.scalar(arrays[label+'_rate'])
        grad=np.stack((grid.r*(grid.dNx@arrays[label]@grid.Nz.T),np.zeros_like(s),grid.Nx@arrays[label]@grid.dNz.T),axis=-1)
        if spec.get('constant_entropy',False):
            ps=f['n']*E*np.exp(2.)*(f['n']/1e20)**(gamma-1);st_exact=np.zeros_like(n)
        else:
            target=ref.thermal(grid.r,grid.z,gamma=gamma);ps=target['pressure']
            gradp=-target['electric']*(E*f['n'])[...,None]
            grad_exact=gradp/ps[...,None]-gamma*f['gradn']/f['n'][...,None]
            st_exact=-np.sum(f['particle']*grad_exact,axis=-1)/f['n']
        pressure=n*E*np.exp(s)*(n/1e20)**(gamma-1);internal=pressure/(gamma-1)
        energy_t=internal*(st+gamma*nt/n);energy_exact=ps/(gamma-1)*(st_exact+gamma*f['nt']/f['n'])
        temperature_t=pressure/(n*E)*(st+(gamma-1)*nt/n)
        temp_exact=ps/(f['n']*E)*(st_exact+(gamma-1)*f['nt']/f['n'])
        metrics[label]=dict(entropy_rate=grid.metric(st,st_exact),pressure=grid.metric(pressure,ps),
                            internal_rate=grid.metric(energy_t,energy_exact),temperature_rate=grid.metric(temperature_t,temp_exact))
        local[label]=grid.cells(energy_t)
        entropy_terms=(n*st,nt*s);weak_terms=(nt*s,np.sum(particle*grad,axis=-1))
        balances[label]=dict(entropy_backward=float(abs(np.sum(grid.weight*(entropy_terms[0]+entropy_terms[1])))/np.sum(grid.weight*(abs(entropy_terms[0])+abs(entropy_terms[1])))),
            weak_continuity_backward=float(abs(np.sum(grid.weight*(weak_terms[0]-weak_terms[1])))/np.sum(grid.weight*(abs(weak_terms[0])+abs(weak_terms[1])))),
            energy_rate=float(np.sum(grid.weight*energy_t)))
    return dict(metrics=metrics,balances=balances),local


def transport_norms(mesh,ref,spec,arrays,order=14):
    grid=PhysicalScalarSampler(mesh,spec.get('degree',5),order);records=[]
    for index,time in enumerate(arrays['times']):
        n=grid.density(arrays['density_snapshots'][index]);exact=ref.fields(grid.r,grid.z,float(time))
        row=dict(time=float(time),density=grid.metric(n,exact['n']),particles=float(np.sum(grid.weight*n)))
        for label,gamma in (('electron',ref.gamma_e),('ion',ref.gamma_i)):
            s=grid.scalar(arrays[label+'_snapshots'][index]);temp=np.exp(s)*(n/1e20)**(gamma-1)
            pressure=n*E*temp;internal=pressure/(gamma-1)
            background=(3+.1*(grid.r/ref.radius)**2) if label=='electron' else (2.5-.15*(grid.r/ref.radius)**2)
            row[label]=dict(entropy_wave=grid.metric(s-background,exact[label]-background),
                temperature=grid.metric(temp,exact[label+'_temperature']),pressure=grid.metric(pressure,exact[label+'_pressure']),
                entropy=float(np.sum(grid.weight*n*s)),internal=float(np.sum(grid.weight*internal)),
                exact_internal=float(np.sum(grid.weight*exact[label+'_internal'])))
        records.append(row)
    return records
