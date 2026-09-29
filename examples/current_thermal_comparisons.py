"""Physical inter-trajectory differences on common or nested meshes."""
import numpy as np


def fields(mesh,case,archive,prefix,r,z):
    r,z=np.broadcast_arrays(r,z)
    i=np.clip(np.searchsorted(mesh.rf,r,side='right')-1,0,mesh.nr-1)
    location=(z+mesh.half_length)/mesh.dz;j=np.floor(location).astype(int)%mesh.nz;tz=location-np.floor(location)
    tx=(r*r-mesh.rf[i]**2)/(mesh.rf[i+1]**2-mesh.rf[i]**2);nextj=(j+1)%mesh.nz
    def scalar(v):return (1-tx)*((1-tz)*v[i,j]+tz*v[i,nextj])+tx*((1-tz)*v[i+1,j]+tz*v[i+1,nextj])
    p=archive[prefix+'_psi'];g=archive[prefix+'_g'];c=archive[prefix+'_current']
    br=-((1-tx)*(p[i,j+1]-p[i,j])+tx*(p[i+1,j+1]-p[i+1,j]))/(r*mesh.dz)
    bz=((1-tz)*(p[i+1,j]-p[i,j])+tz*(p[i+1,j+1]-p[i,j+1]))/mesh.radial_measure[i]
    count=mesh.nr*mesh.nz;theta=(mesh.nr-1)*mesh.nz
    cr=c[:count].reshape(mesh.nr,mesh.nz);ct=np.pad(c[count:count+theta].reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
    cz=c[count+theta:].reshape(mesh.nr+1,mesh.nz)
    jr=r*((1-tz)*cr[i,j]+tz*cr[i,nextj]);jt=scalar(mesh.rf[:,None]*ct)/r;jz=(1-tx)*cz[i,j]+tx*cz[i+1,j]
    result=dict(magnetic=np.stack((br,r*g[i,j],bz),axis=-1),current=np.stack((jr,jt,jz),axis=-1))
    n=archive['number_density']
    for label in ('electron','ion'):
        gamma=case['gamma_e' if label=='electron' else 'gamma_i']
        result[label+'_pressure']=scalar((n/case['density_scale'])**(gamma-1)*archive[prefix+'_'+label+'_entropy'])[...,None]
    return result


def difference(coarse_case,coarse_mesh,coarse,fine_case,fine_mesh,fine,order=6):
    if fine_mesh.nr%coarse_mesh.nr or fine_mesh.nz%coarse_mesh.nz:
        raise ValueError('Comparison requires coincident or nested meshes')
    nodes,weights=np.polynomial.legendre.leggauss(order);errors={};changes={};scales={}
    for a,wa in zip(nodes,weights):
        r=fine_mesh.r[:,None]+.5*fine_mesh.dr*a
        for b,wb in zip(nodes,weights):
            z=fine_mesh.z[None]+.5*fine_mesh.dz*b
            ca=fields(coarse_mesh,coarse_case,coarse,'final',r,z)
            fi=fields(fine_mesh,fine_case,fine,'final',r,z);initial=fields(fine_mesh,fine_case,fine,'initial',r,z)
            weight=.5*np.pi*wa*wb*r*fine_mesh.dr*fine_mesh.dz
            for key in ca:
                errors[key]=errors.get(key,0.)+np.sum(weight*np.sum((ca[key]-fi[key])**2,axis=-1))
                changes[key]=changes.get(key,0.)+np.sum(weight*np.sum((fi[key]-initial[key])**2,axis=-1))
                scales[key]=scales.get(key,0.)+np.sum(weight*np.sum(initial[key]**2,axis=-1))
    return {key:dict(rms_difference=float(np.sqrt(errors[key]/np.sum(fine_mesh.volumes))),
        relative_to_fine_change=float(np.sqrt(errors[key]/changes[key])) if changes[key]>1e-28*scales[key] else None,
        relative_to_initial_field=float(np.sqrt(errors[key]/scales[key])) if scales[key] else None) for key in errors}


def comparisons(records,archives,meshes):
    def pair(a,b):
        assert records[a]['completed'] and records[b]['completed']
        np.testing.assert_allclose(records[a]['final']['time_s'],records[b]['final']['time_s'],rtol=1e-14)
        return difference(records[a]['configuration'],meshes[a],archives[a],records[b]['configuration'],meshes[b],archives[b])
    result={}
    for family,counts,base in (('original',(20,40,80,160),'original_16'),('wave',(10,20,40,80),'wave_h1_8')):
        ref=family+'_time_320';values={}
        for n in counts:
            name=base if n==80 else family+'_time_'+str(n)
            values[str(n)]=pair(name,ref)
        ratios={key:[values[str(a)][key]['rms_difference']/values[str(b)][key]['rms_difference']
            if values[str(b)][key]['rms_difference'] else None for a,b in zip(counts[:-1],counts[1:])]
            for key in values[str(counts[0])]}
        result[family+'_time']=dict(reference_steps=320,errors=values,ratios=ratios)
    result['original_mesh']={f'{n}_to_{2*n}':pair(f'original_{n}',f'original_{2*n}') for n in (8,16,32,64)}
    result['controls']={name:pair(name,base) for name,base in (
        ('original_quadrature','original_16'),('original_scale','original_16'),('original_tolerance','original_16'),
        ('original_fine_half_time','original_128'),('wave_quadrature','wave_h1_16'))}
    return result
