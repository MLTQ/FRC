"""Presentation sampling of archived axisymmetric fields; never advances physics."""
import base64
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from plasma_bench.moving_geometry import MovingDomain


class SavedField:
    """Locate deformed quads with triangles, then evaluate their original Q1 basis."""
    def __init__(self, domain, arrays, prefix):
        self.d=domain
        self.xy=domain.coordinates(arrays[prefix+'_coordinates'])
        self.a=arrays[prefix+'_induced_a'][:domain.node_count]+float(arrays[prefix+'_bias'])/2
        self.q=arrays[prefix+'_current_over_radius']
        triangles=[];owners=[]
        for i,(a,b,c,d) in enumerate(domain.cells):
            for ids in ((a,b,c),(a,c,d)):
                if len(set(ids))==3:
                    triangles.append(ids);owners.append(i)
        self.owners=np.array(owners)
        self.finder=mtri.Triangulation(self.xy[:,0],self.xy[:,1],triangles).get_trifinder()

    def sample(self, points):
        points=np.asarray(points);tri=self.finder(points[:,0],points[:,1])
        valid=tri>=0;out=np.zeros((len(points),5))
        if not valid.any():return out,valid,np.zeros(0,dtype=bool)
        cell=self.owners[tri[valid]];ids=self.d.cells[cell];x=self.xy[ids]
        target=points[valid];uv=np.full((len(cell),2),.5)
        for _ in range(9):
            u,v=uv.T
            shape=np.stack(((1-u)*(1-v),u*(1-v),u*v,(1-u)*v),axis=1)
            deriv=np.stack((np.stack((-(1-v),1-v,v,-v),axis=1),
                            np.stack((-(1-u),-u,u,1-u),axis=1)),axis=-1)
            jac=np.einsum('qnb,qna->qab',deriv,x)
            residual=np.einsum('qn,qna->qa',shape,x)-target
            uv-=np.linalg.solve(jac,residual[...,None])[...,0]
        u,v=uv.T
        shape=np.stack(((1-u)*(1-v),u*(1-v),u*v,(1-u)*v),axis=1)
        deriv=np.stack((np.stack((-(1-v),1-v,v,-v),axis=1),
                        np.stack((-(1-u),-u,u,1-u),axis=1)),axis=-1)
        jac=np.einsum('qnb,qna->qab',deriv,x)
        grad=np.einsum('qnb,qba->qna',deriv,np.linalg.inv(jac))
        aval=np.sum(shape*self.a[ids],axis=1)
        da=np.einsum('qna,qn->qa',grad,self.a[ids]);r=target[:,0]
        br=-r*da[:,1];bz=2*aval+r*da[:,0]
        inside=cell<self.d.fluid_cell_count
        j=np.zeros(len(cell));j[inside]=r[inside]*np.sum(shape[inside]*self.q[ids[inside]],axis=1)
        out[valid]=np.column_stack((j,np.hypot(br,bz),r*r*aval,br,bz))
        assert np.max(np.abs(np.einsum('qn,qna->qa',shape,x)-target))<1e-9
        return out,valid,inside


def build(source=Path('outputs/surface_verification/dynamic/source_g1_32.npz')):
    record=json.loads(source.with_suffix('.json').read_text())
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    assert record['completed'] and digest==record['fields']['sha256']
    c=record['configuration']
    d=MovingDomain(c['radial_cells'],c['angular_cells'],radius=c['radius_m'],
        matching_radius=c['matching_radius'],vacuum_layers=c['vacuum_layers'],
        modes=c['modes'],quadrature=c['quadrature'],radial_grading=c['radial_grading'])
    nr,nz=224,352;rmax,zmax=1.42,1.64
    r=(np.arange(nr)+.5)*rmax/nr;z=((np.arange(nz)+.5)/nz*2-1)*zmax
    rr,zz=np.meshgrid(r,z);points=np.column_stack((rr.ravel(),zz.ravel()))
    frames=[];samples=[];samplers=[];supports=[]
    with np.load(source,allow_pickle=False) as arrays:
        for item in record['fields']['times']:
            prefix=item['prefix'];sampler=SavedField(d,arrays,prefix)
            values,valid,_=sampler.sample(points);values=values.reshape(nz,nr,5)
            sample=values[:,:,:4].copy();sample[:,:,0]=np.abs(sample[:,:,0])
            samples.append(sample);samplers.append(sampler);supports.append(valid.reshape(nz,nr))
            exact_q=arrays[prefix+'_current_quadrature']
            qpoints=arrays[prefix+'_quadrature_coordinates'].reshape(-1,2)
            replay=sampler.sample(qpoints)[0][:,0].reshape(exact_q.shape)
            relative=np.linalg.norm(replay-exact_q)/np.linalg.norm(exact_q)
            assert relative<1e-10
            frames.append(dict(time_ms=item['time_s']*1000,
                bias_mT=float(arrays[prefix+'_bias'])*1000,
                max_current_A_m2=float(np.max(np.abs(exact_q))),
                current_replay_relative=float(relative)))
    jmax=max(s[:,:,0].max() for s in samples);bmax=max(s[:,:,1].max() for s in samples)
    psimin=min(s[:,:,2].min() for s in samples)
    levels=np.r_[np.linspace(psimin*.9,psimin*.08,7),np.linspace(1.5e-5,7.0e-4,10)]
    for frame,sample,sampler,support in zip(frames,samples,samplers,supports):
        fig,ax=plt.subplots();cs=ax.contour(z,r,np.ma.array(sample[:,:,2],mask=~support).T,levels=levels)
        paths=[];err=[]
        for level,segments in zip(levels,cs.allsegs):
            for seg in segments:
                if len(seg)<12:continue
                # At most 200 points per line for a responsive exported viewer.
                indices=np.unique(np.r_[np.arange(0,len(seg),max(1,len(seg)//200)),len(seg)-1])
                line=seg[indices];vals,valid,_=sampler.sample(line[:,::-1])
                if not valid.all():continue
                err.append(float(np.max(np.abs(vals[:,2]-level))/max(abs(psimin),1e-20)))
                paths.append(dict(points=np.round(line,6).tolist(),b=np.round(vals[:,1]*1000,6).tolist(),closed=bool(level<0)))
        plt.close(fig)
        frame['lines']=paths;frame['contour_error_relative_to_peak_flux']=max(err)
        frame['texture']=base64.b64encode(sample.astype('<f4').tobytes()).decode()
        rings=[]
        for rz in [(rad,axial) for axial in (-.48,-.24,0,.24,.48) for rad in (.35,.57,.78,.93)]:
            val,valid,_=sampler.sample(np.array([rz]))
            if valid[0] and abs(val[0,0])>.12*jmax:
                rings.append([*rz,float(abs(val[0,0])/jmax)])
        frame['rings']=rings
    return dict(case=c['name'],source=str(source),source_sha256=digest,
        calculation_source_sha256=record['provenance']['source_sha256'],configuration=c,
        texture_size=[nr,nz],extent=[rmax,zmax],jmax=float(jmax),bmax_mT=float(bmax*1000),frames=frames,
        interpretation='Initialized axisymmetric resistive FRC test; no formation or stability claim.',
        rendering='Exact archived Q1 a and weak J/r evaluated on a display grid; |J| zero outside material. Axisymmetric revolution; poloidal flux contours copied around the axis. Fixed color scales across times. Display interpolation is not a finer simulation.')


if __name__=='__main__':
    output=Path('outputs/frc_showcase');output.mkdir(exist_ok=True)
    data=build();(output/'data.json').write_text(json.dumps(data,separators=(',',':')))
    audit={k:v for k,v in data.items() if k!='frames'}
    audit['frames']=[{k:v for k,v in frame.items() if k not in ('texture','lines','rings')} for frame in data['frames']]
    (output/'data_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit['frames'],indent=2))
